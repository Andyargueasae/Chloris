// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';
import react from '@astrojs/react';

export default defineConfig({
  site: 'https://andyargueasae.github.io',
  base: '/Chloris/',
  integrations: [
    starlight({
      title: 'Chloris',
      disable404Route: false,
      sidebar: [
        { label: 'Home', slug: 'index' },
        { label: 'Description', slug: 'description' },
        { label: 'Species', slug: 'species' },
        { label: 'Phylogeny', slug: 'phylogeny' },
        {
          label: 'Proteins',
          items: [{ autogenerate: { directory: 'proteins' } }],
        },
        { label: 'Download', slug: 'download' },
        { label: 'Credits', slug: 'credits' },
      ],
    }),
    react(),
  ],
});
