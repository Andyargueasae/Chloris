// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// https://astro.build/config
export default defineConfig({
	security: {
		contentSecurityPolicy: {
		'default-src': ["'self'"],
		'script-src':  ["'self'", 'https://unpkg.com', 'https://code.jquery.com', 'https://cdn.jsdelivr.net'],
		'style-src':   ["'self'", 'https://unpkg.com', 'https://cdn.jsdelivr.net', "'unsafe-inline'"],
		'img-src':     ["'self'", 'data:', 'blob:'],
		'connect-src': ["'self'", 'https://unpkg.com', 'https://code.jquery.com', 'https://cdn.jsdelivr.net'],
		'worker-src':  ["'self'", 'blob:'],
		}
},
	integrations: [
		starlight({
			title: 'A2K',
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/Andyargueasae/A2K' },
				{ icon: 'download', label: 'Download', href: '/download' },],
			sidebar: [
				// {
				// 	label: 'Species',
				// 	autogenerate: { directory: 'species' },
				// },
				// {
				// 	label: 'Genes',
				// 	autogenerate: { directory: 'genes' },
				// },
				{ slug: 'description' },
				{ slug: 'species' },
				{ slug: 'genes' },
				{ slug: 'download' },
				{ slug: 'credits' },
			],
		}),
	],
});

