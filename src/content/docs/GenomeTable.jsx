import { useMemo, useState } from 'react'
import raw from '/src/content/galah_cluster_reps.tsv?raw'

const rows = raw.trim().split('\n').map(line => line.split('\t'))
const headers = rows[0].map(h => h.trim())

const data = rows.slice(1).map(row =>
  Object.fromEntries(
    headers.map((h, i) => {
      let val = row[i]?.trim() ?? ""

      if (h === "in_tree") val = val === "True" || val === "true"
      if (h === "num_in_tree") val = Number(val)

      return [h, val]
    })
  )
)

export function GenomeTable() {
  const [query, setQuery] = useState("")
  const [treeFilter, setTreeFilter] = useState("all")
  const [page, setPage] = useState(0)

  const pageSize = 25

  const filtered = useMemo(() => {
    return data.filter(row => {
      const q = query.toLowerCase()

      const matchesQuery =
        row.rep_genome?.toLowerCase().includes(q) ||
        row.genome?.toLowerCase().includes(q) ||
        row.genomes_in_tree?.toLowerCase().includes(q)

      const matchesTree =
        treeFilter === "all" ||
        (treeFilter === "true" && row.in_tree) ||
        (treeFilter === "false" && !row.in_tree)

      return matchesQuery && matchesTree
    })
  }, [query, treeFilter])

  const totalPages = Math.ceil(filtered.length / pageSize)
  const current = filtered.slice(page * pageSize, (page + 1) * pageSize)

  return (
    <div>
      <input
        placeholder="Search genome..."
        value={query}
        onChange={e => {
          setQuery(e.target.value)
          setPage(0)
        }}
      />

      <select
        value={treeFilter}
        onChange={e => {
          setTreeFilter(e.target.value)
          setPage(0)
        }}
      >
        <option value="all">All genomes</option>
        <option value="true">In tree only</option>
        <option value="false">Not in tree</option>
      </select>

      <p>Showing {current.length} of {filtered.length} genomes</p>

      <div style={{ maxHeight: 600, overflow: "auto" }}>
        <table>
          <thead>
            <tr>
              <th>Representative genome</th>
              <th>Genome</th>
              <th>In tree</th>
              <th>Number in tree</th>
              <th>Genomes in tree</th>
              <th>GenBank</th>
            </tr>
          </thead>

          <tbody>
            {current.map((row, i) => (
              <tr key={`${row.rep_genome}-${i}`}>
                <td>{row.rep_genome}</td>
                <td>{row.genome}</td>
                <td>{row.in_tree ? "True" : "False"}</td>
                <td>{row.num_in_tree}</td>
                <td>{row.genomes_in_tree || "—"}</td>
                <td>
                  <a href={row.GenBank_URL} target="_blank" rel="noopener noreferrer">
                    GenBank
                  </a>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>

      <button onClick={() => setPage(p => Math.max(p - 1, 0))} disabled={page === 0}>
        Previous
      </button>

      <span> Page {page + 1} of {totalPages || 1} </span>

      <button
        onClick={() => setPage(p => Math.min(p + 1, totalPages - 1))}
        disabled={page >= totalPages - 1}
      >
        Next
      </button>
    </div>
  )
}