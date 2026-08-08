import React from 'react';
import App from './App';
import { createRoot } from 'react-dom/client';

if (import.meta.env.PROD) {
    if ('serviceWorker' in navigator) {
        navigator.serviceWorker.register(
            import.meta.env.BASE_URL + 'sw.js',
            { type: 'module' }
        )
    }
}

createRoot(document.getElementById('app')!).render(<App />);