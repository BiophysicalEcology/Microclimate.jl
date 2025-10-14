import { defineConfig } from 'vitepress'

export default defineConfig({
  title: "Microclimate.jl",
  description: "A Julia package for simulating microclimates.",
  base: '/Microclimate.jl/',
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: 'Home', link: '/' },
      { text: 'API', link: '/api' }
    ],

    sidebar: [
      {
        text: 'Examples',
        items: [
          { text: 'Item A', link: '/item-a' },
          { text: 'Item B', link: '/item-b' }
        ]
      }
    ],

    socialLinks: [
      { icon: 'github', link: 'https://github.com/BiophysicalEcology/Microclimate.jl' }
    ]
  }
})
