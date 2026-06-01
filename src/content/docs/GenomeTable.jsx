import { useEffect, useMemo, useState } from 'react'
import { AgGridReact } from 'ag-grid-react'
import { AllCommunityModule, ModuleRegistry } from 'ag-grid-community'
import 'ag-grid-community/styles/ag-grid.css'
import 'ag-grid-community/styles/ag-theme-alpine.css'

ModuleRegistry.registerModules([AllCommunityModule])

export function GenomeTable() {
  const [data, setData] = useState([])
  const [isLoading, setIsLoading] = useState(true)
  const [error, setError] = useState(null)
  const [quickFilterText, setQuickFilterText] = useState('')

  const defaultColDef = useMemo(() => {
    return {
      sortable: true,
      filter: true,
      floatingFilter: true,
      resizable: true,
      editable: false,
      minWidth: 140,
      flex: 1,
    }
  }, [])

  const columnDefs = useMemo(() => {
    return [
      {
        field: 'rep_genome',
        headerName: 'Representative Genome',
        pinned: 'left',
        minWidth: 220,
        filter: 'agTextColumnFilter',
      },
      {
        field: 'species',
        headerName: 'Species',
        minWidth: 220,
        filter: 'agTextColumnFilter',
      },
      {
        field: 'members',
        headerName: 'Members',
        minWidth: 360,
        flex: 2,
        filter: 'agTextColumnFilter',
        cellStyle: {
          whiteSpace: 'normal',
          lineHeight: '1.4',
        },
        autoHeight: true,
      },
      {
        field: 'genomes_in_tree',
        headerName: 'Genomes in Tree',
        minWidth: 200,
        filter: 'agTextColumnFilter',
        cellRenderer: params => {
          const value = params.value
          if (!value) {
            return <span style={{ color: '#999' }}>-</span>
          }
          return (
            <a
              href={`https://www.ncbi.nlm.nih.gov/nuccore/${value}`}
              target="_blank"
              rel="noopener noreferrer"
              style={{ color: '#0066cc', textDecoration: 'none', fontWeight: 500 }}
              onMouseEnter={e => {
                e.currentTarget.style.textDecoration = 'underline'
              }}
              onMouseLeave={e => {
                e.currentTarget.style.textDecoration = 'none'
              }}
            >
              {value}
            </a>
          )
        },
      },
    ]
  }, [])

  const sideBar = useMemo(() => {
    return {
      toolPanels: ['columns', 'filters'],
      defaultToolPanel: 'columns',
    }
  }, [])

  useEffect(() => {
    const loadData = async () => {
      try {
        const tsvUrl = `${import.meta.env.BASE_URL}downloads/galah_cluster_reps_filtered.tsv`
        const response = await fetch(tsvUrl)

        if (!response.ok) {
          throw new Error(`HTTP error! status: ${response.status}`)
        }

        const raw = await response.text()

        const rows = raw.trim().split('\n').map(line => line.split('\t'))
        const headers = rows[0].map(h => h.trim())

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

        setData(parsedData)
        setError(null)
      } catch (err) {
        setError(err.message)
      } finally {
        setIsLoading(false)
      }
    }

    loadData()
  }, [])

  if (isLoading) {
    return <p>Loading genome data...</p>
  }

  if (error) {
    return (
      <div style={{ padding: '1rem', backgroundColor: '#ffe0e0', color: '#c00', borderRadius: '4px' }}>
        <p><strong>Error:</strong> {error}</p>
      </div>
    )
  }

  if (data.length === 0) {
    return <p>No genome data available</p>
  }

  return (
    <div style={{ marginTop: '2rem', width: '100%' }}>
      <p style={{ marginBottom: '1rem', color: '#333', fontWeight: '500' }}>
        Showing {data.length} genomes
      </p>

      <div style={{ display: 'flex', gap: '0.75rem', marginBottom: '0.75rem', alignItems: 'center' }}>
        <label htmlFor="genome-quick-filter" style={{ fontWeight: 500, color: '#333' }}>
          Quick filter:
        </label>
        <input
          id="genome-quick-filter"
          type="text"
          value={quickFilterText}
          onChange={e => setQuickFilterText(e.target.value)}
          placeholder="Type to filter all columns"
          style={{
            padding: '0.45rem 0.6rem',
            border: '1px solid #c8c8c8',
            borderRadius: '4px',
            minWidth: '300px',
            maxWidth: '100%',
          }}
        />
      </div>

      <div
        className="ag-theme-alpine"
        style={{
          width: '100%',
          height: '75vh',
          minHeight: '650px',
          border: '1px solid #d5d5d5',
          borderRadius: '6px',
          overflow: 'hidden',
        }}
      >
        <AgGridReact
          rowData={data}
          columnDefs={columnDefs}
          defaultColDef={defaultColDef}
          quickFilterText={quickFilterText}
          sideBar={sideBar}
          animateRows={true}
          enableCellTextSelection={true}
          pagination={true}
          paginationPageSize={25}
          paginationPageSizeSelector={[10, 25, 50, 100]}
          rowSelection={{ mode: 'multiRow' }}
          suppressRowClickSelection={false}
        />
      </div>
    </div>
  )
}