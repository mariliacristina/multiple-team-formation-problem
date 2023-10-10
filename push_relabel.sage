# Calcula o fluxo maximo de s para t, chamando a operacao de discharge para cada vertice em L
def relabel_to_front(G, s, t):

	# Dicionario para armazenar atributos dos vertices: altura, excesso e lista de vizinhos, nessa ordem
	height = dict()
	excess = dict()
	neighbors = dict()

	# Dicionario para armazenar a quantidade de fluxo em cada aresta 
	flow = dict()

	# Dicionario que representa a rede residual, contendo a capacidade residual de cada aresta
	residual = dict()
	
	# Preenche a lista de vizinhos para cada vertice de G
	def fill_neighbors():
		for u in G.vertices():
			neighbors[u] = G.neighbors(u)
	
	fill_neighbors()
	
	# Cria um pre-fluxo inicial na rede de fluxo
	def initialize_preflow():
		# Para cada vertice em G, inicializa a altura e o excesso com 0, respectivamente
		for v in G.vertices():
			height[v] = 0
			excess[v] = 0
		
		# Para cada aresta em G, inicializa o fluxo com 0
		for (u, v, c) in G.edges():
			flow[(u, v)] = 0
			
			residual[(u, v)] = c
			residual[(v, u)] = 0
		
		# Altura de s = numero de vertices de G
		height[s] = G.num_verts()
		
		# Para cada v adjacente a s
		for (_, v, c) in G.outgoing_edges(s):
			flow[(s, v)] = c
			
			residual[(s, v)] = 0
			residual[(v, s)] = c
			
			excess[v] = c
			excess[s] -= c

	initialize_preflow()
	
	# Empurra fluxo de u para v
	def push(u, v):
		delta = min(excess[u], residual[(u, v)])
		
		if G.has_edge(u, v):
			flow[(u, v)] += delta
			
			residual[(u, v)] -= delta
			residual[(v, u)] += delta
		
		else:
			flow[(v, u)] -= delta
			
			residual[(v, u)] += delta
			residual[(u, v)] -= delta
		
		excess[u] -= delta
		excess[v] += delta


	# Remarca o vertice u
	def relabel(u):
		height[u] = 1 + min(height[v] for v in G.neighbors(u) if residual[(u, v)] > 0)


	# Descarrega (elimina) todo o excesso de fluxo do vertice u
	def discharge(u):
		i = 0
		
		while excess[u] > 0:
			if i < len(neighbors[u]):
				v = neighbors[u][i]
				
				if (residual[(u, v)] > 0) and (height[u] == height[v] + 1):
					push(u, v)
				
				else:
					i += 1
			
			else:
				relabel(u)
				i = 0
	
	# L = G.V - {s,t}
	# G.vertices() esta ordenado pelo nome dos vertices em ordem alfabetica. Como so temos os vertices 'Hi', 'Rij', 's' e 't', os dois ultimos sempre serao 's' e 't', logo o comando abaixo seleciona exatamente G.V - {s,t}.
	L = G.vertices()[:-2]
	
	i = 0
	
	while i < len(L): 
		u = L[i]
		
		old_height = height[u]
		
		discharge(u)
		
		# Mover u para a frente da lista L, e voltar o algoritmo para o segundo vertice de L (i = 1) 
		if height[u] > old_height:
			L.remove(u)
			L = [u] + L
			i = 1
		
		else:
			i += 1


	return flow
	

