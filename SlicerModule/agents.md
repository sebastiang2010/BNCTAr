# Antigravity Instructions

Este proyecto está conectado a:

<https://github.com/sebastiang2010/BNCTAr.git>

## Reglas

* usar commits descriptivos
* realizar push automático a origin/main
* revisar sintaxis Python antes de commit
* priorizar estabilidad sobre optimización

## Nunca subir

* .env
* API keys
* logs
* **pycache**
* credenciales Binance
* archivos temporales

## Formato de commits

fix: corrige sincronización websocket
feat: agrega trailing stop
refactor: mejora reconciliación REST/ws

## Trading Bot

Este proyecto es un bot HFT/scalper para Binance Futures.

Tener especial cuidado con:

* fills tardíos
* sincronización websocket/rest
* manejo de posiciones
* cancelaciones de órdenes
* race conditions
