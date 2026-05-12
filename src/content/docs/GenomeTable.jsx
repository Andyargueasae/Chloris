import { useEffect, useState } from 'react'

export function GenomeTable() {
  const [data, setData] = useState([])
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)

  useEffect(() => {
    const loadData = async () => {
      try {
        console.log('🔍 GenomeTable component mounted')
        const baseUrl = import.meta.env.BASE_URL || '/'
        const tsvUrl = `${baseUrl}galah_cluster_reps_filtered.tsv`
        const response = await fetch(tsvUrl)
        console.log('Fetch response:', response.status)
        
        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`)
        }
        
        const raw = await response.text()
        console.log('Raw data length:', raw.length)
        
        const rows = raw.trim().split('\n').map(line => line.split('\t'))
        const headers = rows[0].map(h => h.trim())
        console.log('Headers:', headers)

        const parsedData = rows.slice(1).map(row =>
          Object.fromEntries(
            headers.map((h, i) => {
              let val = row[i]?.trim() ?? ""
              if (h === "in_tree") val = val === "True" || val === "true"
              if (h === "num_in_tree") val = Number(val)
              return [h, val]
            })
          )
        )
        
        console.log('✅ Loaded', parsedData.length, 'genomes')
        setData(parsedData)
        setError(null)
      } catch (err) {
        console.error('❌ Error:', err.message)
        setError(err.message)
      } finally {
        setIsLoading(false)
      }
    }

    loadData()
  }, [])

  if (isLoading) {
    return <p>⏳ Loading genome data...</p>
  }

  if (error) {
    return (
      <div style={{ padding: '1rem', backgroundColor: '#ffe0e0', color: '#c00', borderRadius: '4px' }}>
        <p><strong>Error:</strong> {error}</p>
      </div>
    )
  }

  if (data.length === 0) {
    return <p>⚠️ No genome data available</p>
  }

  return (
    <div style={{ marginTop: '2rem', width: '100%' }}>
      <p style={{ marginBottom: '1rem', color: '#333', fontWeight: '500' }}>
        Showing all {data.length} genomes
      </p>

      <div style={{ overflowX: 'auto', marginBottom: '1.5rem', border: '1px solid #ccc', borderRadius: '4px', maxHeight: '1200px', overflowY: 'auto' }}>
        <table style={{
          width: '100%',
          borderCollapse: 'collapse',
          fontSize: '14px',
          backgroundColor: '#fff',
        }}>
          <thead>
            <tr style={{ backgroundColor: '#0066cc', color: '#fff', fontWeight: '600', position: 'sticky', top: 0 }}>
              <th style={{ padding: '12px', textAlign: 'left', borderRight: '1px solid #ddd' }}>Representative Genome</th>
              <th style={{ padding: '12px', textAlign: 'left', borderRight: '1px solid #ddd' }}>Members</th>
              <th style={{ padding: '12px', textAlign: 'left', borderRight: '1px solid #ddd' }}>Species</th>
              <th style={{ padding: '12px', textAlign: 'left' }}>Genomes in Tree</th>
            </tr>
          </thead>
          <tbody>
            {data.map((row, i) => (
              <tr key={i} style={{ borderBottom: '1px solid #ddd', backgroundColor: i % 2 === 0 ? '#f9f9f9' : '#fff' }}>
                <td style={{ padding: '12px', borderRight: '1px solid #ddd', color: '#333', fontWeight: '500' }}>
                  {row.rep_genome}
                </td>
                <td style={{ padding: '12px', borderRight: '1px solid #ddd', color: '#555', fontSize: '13px', maxWidth: '300px', wordBreak: 'break-word' }}>
                  {row.members}
                </td>
                <td style={{ padding: '12px', borderRight: '1px solid #ddd', color: '#333' }}>
                  {row.species || <span style={{ color: '#999' }}>—</span>}
                </td>
                <td style={{ padding: '12px', color: '#333' }}>
                  {row.genomes_in_tree ? (
                    <a 
                      href={`https://www.ncbi.nlm.nih.gov/nuccore/${row.genomes_in_tree}`}
                      target="_blank" 
                      rel="noopener noreferrer"
                      style={{ 
                        color: '#0066cc', 
                        textDecoration: 'none',
                        fontWeight: '500',
                      }}
                      onMouseEnter={(e) => e.target.style.textDecoration = 'underline'}
                      onMouseLeave={(e) => e.target.style.textDecoration = 'none'}
                    >
                      {row.genomes_in_tree}
                    </a>
                  ) : (
                    <span style={{ color: '#999' }}>—</span>
                  )}
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  )
}