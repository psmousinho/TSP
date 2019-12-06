# Travelling Salesman Problem:
Resolução do TSP de forma exata utilizando o seguint modelo:

<img src="https://i.imgur.com/IyLIeUq.png" title="Model" />


## Execução:
`python sintatico.py [-l -m -rd] {arquivo de entrada}`

Onde:

-l: Salva o log da execução do arquivo results/log

-m: Reproduz o modelo utilizado no cplex no arquivo results/modelo.lp

-rd: Salva a saida do programa de forma tabular em markdown no arquivo README.md

O arquivo de entrada é qualquer instancia paro o problemas tsp da tsplib: http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsplib.html.

## Dependências:
- Python
- Cplex

## Resultados:
|Instância|Status da Solucao|Custo da Solucao|Duracao(seg)|
|---------|-----------------|----------------|------------|
|../instances/burma14.tsp|MIP_optimal|3323.000000000014|1.9648902416229248|
|../instances/ulysses16.tsp|MIP_optimal|6858.999999999987|10.700449705123901|
|../instances/gr17.tsp|MIP_optimal|2084.999999999999|5.252423524856567|
|../instances/gr21.tsp|MIP_optimal|2706.999999999992|1.9691967964172363|
|../instances/ulysses22.tsp|MIP_optimal|7012.999999993599|46.9362850189209|
|../instances/gr24.tsp|MIP_optimal|1272.0000000000418|125.89890050888062|
|../instances/fri26.tsp|MIP_optimal|937.0000000000975|23.447660207748413|
|../instances/bayg29.tsp|MIP_optimal|1609.9999999999927|529.0106751918793|
|../instances/bays29.tsp|MIP_optimal|2020.0000000000557|307.5237355232239|
