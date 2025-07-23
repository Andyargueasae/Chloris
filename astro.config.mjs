// @ts-check
import { defineConfig } from 'astro/config';
import starlight from '@astrojs/starlight';

// https://astro.build/config
export default defineConfig({
	integrations: [
		starlight({
			title: 'A2K',
			social: [
				{ icon: 'github', label: 'GitHub', href: 'https://github.com/withastro/starlight' },
				{ icon: 'discord', label: 'Discord', href: 'https://astro.build/chat' },],
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
