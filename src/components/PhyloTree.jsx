import { useEffect, useRef, useState } from 'react';
import * as d3 from 'd3';

const RANKS = [
  ['d__', 'Domain'],
  ['p__', 'Phylum'],
  ['c__', 'Class'],
  ['o__', 'Order'],
  ['f__', 'Family'],
  ['g__', 'Genus'],
  ['s__', 'Species'],
];

function parseTaxon(name = '') {
  const parts = name.split('_');
  const result = {};

  let currentRank = null;
  let buffer = [];

  for (const part of parts) {
    const rank = RANKS.find(([prefix]) => part.startsWith(prefix));

    if (rank) {
      if (currentRank) result[currentRank] = buffer.join('_');
      currentRank = rank[0];
      buffer = [part.replace(rank[0], '')];
    } else if (currentRank) {
      buffer.push(part);
    }
  }

  if (currentRank) result[currentRank] = buffer.join('_');

  return result;
}

function displayLabel(name = '') {
  const tax = parseTaxon(name);

  return (
    tax['s__'] ||
    tax['g__'] ||
    tax['f__'] ||
    tax['o__'] ||
    tax['c__'] ||
    tax['p__'] ||
    tax['d__'] ||
    name
  );
}

function parseNewick(newick) {
  const ancestors = [];
  const tree = {};
  let current = tree;

  newick.split(/\s*(;|\(|\)|,|:)\s*/).forEach((token, i, tokens) => {
    if (!token) return;

    switch (token) {
      case '(': {
        const child = {};
        current.children = [child];
        ancestors.push(current);
        current = child;
        break;
      }

      case ',': {
        const sibling = {};
        ancestors[ancestors.length - 1].children.push(sibling);
        current = sibling;
        break;
      }

      case ')':
        current = ancestors.pop();
        break;

      case ':':
      case ';':
        break;

      default: {
        const previous = tokens[i - 1];

        if (previous === ':' && !Number.isNaN(Number(token))) {
          current.length = Number(token);
        } else {
          current.name = token;
        }
      }
    }
  });

  return tree;
}

function countLeaves(node) {
  if (!node.children && !node._children) return 1;
  const children = node.children || node._children || [];
  return children.reduce((sum, child) => sum + countLeaves(child), 0);
}

function getCommonRankLabel(node) {
  const leaves = node.leaves ? node.leaves() : [];
  if (!leaves.length) return '';

  for (const [prefix, rankName] of RANKS) {
    const values = new Set(
      leaves
        .map(leaf => parseTaxon(leaf.data.name || '')[prefix])
        .filter(Boolean)
    );

    if (values.size === 1) {
      const value = [...values][0];
      return `${rankName}: ${value}`;
    }
  }

  return `${leaves.length} taxa`;
}

export default function PhyloTree({ treeUrl }) {
  const containerRef = useRef(null);
  const [error, setError] = useState('');

  useEffect(() => {
    async function loadTree() {
      try {
        const res = await fetch(treeUrl);
        if (!res.ok) throw new Error(`Could not fetch tree: ${res.status}`);

        const newick = await res.text();
        const data = parseNewick(newick.trim());

        containerRef.current.innerHTML = '';

        const width = 1800;
        const dx = 14;
        const dy = 160;

        const root = d3.hierarchy(data);

        root.each(d => {
          d._leafCount = countLeaves(d);
          d._commonLabel = getCommonRankLabel(d);
        });

        function collapseByDepth(node, maxDepth = 4) {
          if (node.depth >= maxDepth && node.children) {
            node._children = node.children;
            node.children = null;
          }

          (node.children || node._children || []).forEach(child =>
            collapseByDepth(child, maxDepth)
          );
        }

        collapseByDepth(root, 4);

        const treeLayout = d3.tree().nodeSize([dx, dy]);

        const svg = d3
          .select(containerRef.current)
          .append('svg')
          .attr('width', width)
          .attr('height', 900)
          .style('max-width', '100%')
          .style('height', 'auto')
          .style('font-family', 'system-ui, sans-serif');

        const g = svg.append('g').attr('transform', 'translate(60,60)');

        svg.call(
          d3.zoom().on('zoom', event => {
            g.attr('transform', event.transform);
          })
        );

        function update() {
          treeLayout(root);

          const nodes = root.descendants();
          const links = root.links();

          const minX = d3.min(nodes, d => d.x);
          const maxX = d3.max(nodes, d => d.x);
          const maxY = d3.max(nodes, d => d.y);

          svg
            .attr('height', Math.max(900, maxX - minX + 160))
            .attr('width', Math.max(1200, maxY + 500));

          g.selectAll('*').remove();

          g.append('g')
            .attr('fill', 'none')
            .attr('stroke', '#9ca3af')
            .attr('stroke-width', 1)
            .selectAll('path')
            .data(links)
            .join('path')
            .attr(
              'd',
              d =>
                `M${d.source.y},${d.source.x}
                 H${d.target.y}
                 V${d.target.x}`
            );

          const node = g
            .append('g')
            .selectAll('g')
            .data(nodes)
            .join('g')
            .attr('transform', d => `translate(${d.y},${d.x})`)
            .style('cursor', d =>
              d.children || d._children ? 'pointer' : 'default'
            )
            .on('click', (event, d) => {
              event.stopPropagation();

              if (d.children) {
                d._children = d.children;
                d.children = null;
              } else if (d._children) {
                d.children = d._children;
                d._children = null;
              }

              update();
            });

          node
            .append('circle')
            .attr('r', d => {
              if (d._children) return Math.min(16, 5 + Math.sqrt(d._leafCount));
              return d.children ? 3.5 : 2;
            })
            .attr('fill', d => {
              if (d._children) return '#2563eb';
              if (d.children) return '#64748b';
              return '#334155';
            })
            .attr('stroke', '#ffffff')
            .attr('stroke-width', 1);

          node
            .filter(d => d._children)
            .append('text')
            .attr('dy', '0.35em')
            .attr('text-anchor', 'middle')
            .style('font-size', '9px')
            .style('font-weight', '600')
            .style('fill', '#ffffff')
            .style('pointer-events', 'none')
            .text(d => d._leafCount);

          node
            .filter(d => d._children)
            .append('text')
            .attr('x', 20)
            .attr('dy', '0.31em')
            .style('font-size', '11px')
            .style('font-weight', '600')
            .style('fill', '#1e40af')
            .text(d => d._commonLabel);

          node
            .filter(d => !d.children && !d._children)
            .append('text')
            .attr('dy', '0.31em')
            .attr('x', 7)
            .style('font-size', '10px')
            .style('fill', '#374151')
            .text(d => displayLabel(d.data.name || ''));

          node.append('title').text(d => {
            if (d._children) return `${d._commonLabel}\n${d._leafCount} leaves`;
            return d.data.name || '';
          });
        }

        update();
      } catch (err) {
        console.error(err);
        setError(err.message);
      }
    }

    loadTree();
  }, [treeUrl]);

  return (
    <>
      {error && <pre style={{ color: 'red' }}>{error}</pre>}

      <p style={{ fontSize: '0.9rem', color: '#475569' }}>
        Click blue clades to expand/collapse. Hover labels to see full taxonomy.
      </p>

      <div
        ref={containerRef}
        style={{
          width: '100%',
          height: '900px',
          overflow: 'auto',
          border: '1px solid #d1d5db',
          borderRadius: '8px',
          background: '#ffffff',
        }}
      />
    </>
  );
}