import random
import copy
import math
import time

# Carrega o arquivo que implementa o algoritmo push_relabel
load("/home/marilia/Dropbox/TCC - AGRUPAMENTO/Implementacao/push_relabel.sage")

# Le o arquivo S.txt, que representa a matriz sociometrica S
def read_S(file_name):
	s_file = open(file_name, 'r')
	
	# n = numero de pessoas
	n = int(s_file.readline().strip())
	
	S = []
	
	for line in s_file:
		S.append(line.split())

	for i in range(0,len(S)):
		for j in range(0, len(S[i])):
			S[i][j] = int(S[i][j])

	s_file.close()

	return (S, n)


# Le o arquivo K.txt, que contem as habilidades das pessoas
def read_K(file_name):
	k_file = open(file_name, 'r')
	
	# f = numero de habilidades
	f = int(k_file.readline().strip())
	
	K = []
	
	for line in k_file:
		K.append(line.split())
	
	for i in range(0,len(K)):
		for j in range(0, len(K[i])):
			K[i][j] = int(K[i][j])
	
	k_file.close()

	return (K, f)


# Le o arquivo R.txt, que contem as demandas de habilidades dos projetos
def read_R(file_name):
	r_file = open(file_name, 'r')
	
	# m = numero de projetos
	m = int(r_file.readline().strip())

	R = []
	
	for line in r_file:
		R.append(line.split())
	
	for i in range(0,len(R)):
		for j in range(0, len(R[i])):
			R[i][j] = float(R[i][j])
	
	r_file.close()

	return (R, m)
	
	
# Le o arquivo D.txt, que representa o conjunto D de particionamento de tempos possiveis
def read_D(file_name):
	d_file = open(file_name, 'r')
	
	# t = tamanho de D
	t = int(d_file.readline().strip())

	D = []

	for line in d_file:
		D = line.split()
	
	for i in range(0,len(D)):
			D[i] = float(D[i])
	
	# alfa = diferenca entre os elementos em D
	alfa = 0
	
	if(len(D) > 1):
		alfa = D[1] - D[0]
	else:
		alfa = D[0]

	d_file.close()

	return (D, t, alfa)


# Le o arquivo W.txt, que contem as prioridades de cada projeto
def read_W(file_name):
	w_file = open(file_name, 'r')
	
	# p = numero de projetos
	p = int(w_file.readline().strip())

	W = []

	for line in w_file:
		W = line.split()
	
	for i in range(0,len(W)):
			W[i] = float(W[i])

	w_file.close()

	return (W, p)


# Constroi a rede para obter uma solucao inicial a partir do fluxo maximo       
def construct_network(n, R, alfa, K):
	G = DiGraph()
	
	# Adicionando os vertices fonte (s) e sumidouro (t)
	G.add_vertex('s')
	G.add_vertex('t')
	
	# Adicionando um vertice para cada pessoa
	for i in range(0, n):
		G.add_vertex('H-' + str(i))
		
		# Adicionando uma aresta saindo de s e chegando em cada vertice de uma pessoa, com capacidade 1/alfa
		G.add_edge('s', 'H-' + str(i), 1/alfa)
	
	# arrayP armazena o projeto, arrayS armazena a habilidade requerida pelo projeto
	arrayP = []
	arrayS = []
	
	# Adicionando um vertice para cada demanda de habilidade de cada projeto
	for i in range(0, len(R)):
		for j in range(0, len(R[i])):
			if(R[i][j] != 0.0):
				G.add_vertex('R-' + str(i) + '-' + str(j))
				
				# Adicionando uma aresta saindo de cada vertice de demanda de habilidade e chegando em t com capacidade R[i][j]/alfa
				G.add_edge('R-' + str(i) + '-' + str(j), 't', R[i][j]/alfa)
				
				arrayP.append(i)
				arrayS.append(j)
				
	# Adicionando arestas saindo das pessoas e chegando nas demandas de habilidades
	num_edges_H_R = 0
	for i in range(0, n):
		for j in range(0, len(arrayP)):
			
			# Se a pessoa i possui demanda da habilidade arrayS[j]
			if(K[i][arrayS[j]] != 0):
				G.add_edge('H-' + str(i), 'R-' + str(arrayP[j]) + '-' + str(arrayS[j]), 1/alfa)
				
				num_edges_H_R += 1
				
	return (G, num_edges_H_R)


# Calcula o fluxo maximo de G, com o algoritmo de FF 
def max_flow_FF(G):
	return G.flow('s', 't', value_only=False)


# Verifica se a solucao é viável
def is_feasible(flow, G):
	# A solucao é viável se as arestas que chegam em 't' estão saturadas; 
	# i = (u, v)
	for i in flow:
		u = i[0]
		v = i[1]
		
		if v == 't':
			if flow[i] != G.edge_label(u, v):
				return False
	
	return True


# Constroi uma solucao para o problema, a partir do grafo G
def build_solution(G, m, f, alfa, D):
	# Constroi a solucao a partir do fluxo maximo em G, usando o algoritmo push_relabel de Goldberg
	flow = relabel_to_front(G, 's', 't')
	
	if is_feasible(flow, G) == False:
		return -1
	
	# solution sera da forma: {ProjetoP: {HabilidadeS: [[PessoaH, TempoT]]}}
	# solution sera da forma: {ProjetoP: {HabilidadeS: {TempoT: [Pessoa1, Pessoa2...]}}}
	solution = dict()
	
	# Para cada projeto p, cria-se um subgrupo para cada habilidade s e para cada habilidade s, cria-se um subgrupo para cada
	# tempo t e neste havera um array de pessoas. 
	for p in range(0, m):
		solution[p] = dict()
		
		for s in range(0, f): 
			#solution[p][s] = []
			solution[p][s] = dict()
			
			for t in D:
				solution[p][s][t] = [] # array de pessoas trabalhando no projeto p, com a habilidade s, no tempo t
				
	
	for j in flow:
		u = j[0].split('-')
		v = j[1].split('-')
		
		# Se for uma aresta (Hi, Rjk), o fluxo que passa nela diz a quantidade de tempo alocada da pessoa i no projeto j com a 
		# habilidade k
		if u[0] == 'H' and v[0] == 'R' and flow[j] > 0.0:
			p = int(v[1]) # indice do projeto
			s = int(v[2]) # indice da habilidade
			h = int(u[1]) # indice da pessoa
			
			t = flow[j]*alfa # tempo exercido 
			
			# flow[j] * alfa = fraçao de tempo destinada ao projeto
			solution[p][s][t].append(h) 
			#solution[p][s].append([h, flow[j]*alfa]) 
			#solution[p].append((h, s, flow[j])) 

	
	# TESTE PARA COMPARAR FF COM PUSH_RELABEL 
	flow_golberg = 0
	for i in flow:
		if i[0] == 's':
			flow_golberg += flow[i]
			
	#print("Push-relabel = " + str(flow_golberg))
	#print("FF = " + str(max_flow_FF(G)[0]))
	
	#if(flow_golberg != max_flow_FF(G)[0]):
		#print("Ops! Resultados do FF e do push-relabel não bateram!")
 
	return solution


# Constroi um conjunto inicial de solucoes
def initials_solutions(S, K, R, D, n, f, m, t, alfa):
	# Constroi a rede G
	G, num_edges_H_R = construct_network(n, R, alfa, K)
	
	# Array de solucoes iniciais viaveis
	solutions = []
	
	# Gera a primeira solução
	solut = build_solution(G, m, f, alfa, D)
	
	# Verifica se a solução é viável
	if solut != -1:
		solutions.append(solut)
	
	# Constroi 50 solucoes viaveis 
	for i in range(0,50):
		#Será sorteado um arco para diminuir sua capacidade e gerar outra solucao	
		# Sorteio de um arco (Hi, Rij) aleatorio; 
		# As arestas em G.edges() estão ordenadas da seguinte maneira: (Hi, Rij), (Rij, t), (s, Hi). Então, ao selecionar um
		# número entre 0 e num_edges_H_R, será selecionada exatamente uma aresta (Hi, Rij) 
		edge_index = random.randint(0, num_edges_H_R - 1)
		edge = G.edges()[edge_index]			
		
		# Sorteio de uma capacidade aleatoria
		cap_index = random.randint(0, len(D) - 1)
		cap = D[cap_index]
		
		# Diminui a capacidade do arco em cap/alfa unidades 
		old_cap = edge[2]
		G.set_edge_label(edge[0], edge[1], old_cap - (cap/alfa))
		
		# Constroi uma nova solucao em G, e adiciona ao array solutions
		solut = build_solution(G, m, f, alfa, D)
		
		# Verifica se a solução é viável
		if solut != -1:
			solutions.append(solut)
		
		# O arco volta para sua capacidade original
		G.set_edge_label(edge[0], edge[1], old_cap)	
		
		
	return solutions


# Avalia a qualidade de um projeto
#def project_efficiency(solution, p, S, R, alfa, W):
def project_efficiency(solution, p, S, R):
	e = 0
		
	# Array de pessoas do projeto p, da forma: [[PessoaX, TempoX], [PessoaY, TempoY],...]
	people = []
	
	# Para cada habilidade s
	for s in solution[p]:
		
		# Para cada tempo t
		for t in solution[p][s]:
		
			# Para cada pessoa h	
			for h in solution[p][s][t]:
				people.append([h, t])
	
	# Para cada par de pessoas Hi e Hj 
	for i in people:
		Hi = i[0] # indice da pessoa i
		Xip = i[1] # indice do tempo de i no projeto p
		
		for j in people:
			Hj = j[0] # indice da pessoa j
			Xjp = j[1] # indice do tempo de j no projeto p
			
			#if Hi != Hj:
			e += S[Hi][Hj] * Xip * Xjp
	
	Rsp = 0
	
	# Para cada demanda de habilidade no projeto p
	for r in R[p]:
		Rsp += r
	
	# Calcula a eficiencia do projeto p 
	#efficiency = (((e / Rsp**2) + 1) / 2) 
	efficiency = (((e / Rsp**2) + 1) / 2) 
	
	return efficiency


# Avalia a qualidade de uma solucao
#def global_efficiency(solution, S, R, alfa, W):
def global_efficiency(solution, S, R, m):
	# Eficiencia Global
	efficiency = 0
	
	# Para cada projeto p, acumula a eficiencia dele na eficiencia global
	for p in solution:
		proj_efficiency = project_efficiency(solution, p, S, R)
		
		# efficiency += proj_efficiency * W[p]
		efficiency += proj_efficiency * (1/m)
		
	return efficiency


'''
# Avalia a qualidade de uma solucao
#def global_efficiency(solution, S, R, W):
def global_efficiency(solution, S, R, m):
	# Eficiencia Global
	efficiency = 0
	
	# Para cada projeto p
	for p in solution:
		e = 0
		
		# Array de pessoas do projeto p, da forma: [[PessoaX, TempoX], [PessoaY, TempoY],...]
		people = []
		
		# Para cada habilidade s
		for s in solution[p]:
			
			# Para cada tempo t
			for t in solution[p][s]:
			
				# Para cada pessoa h	
				for h in solution[p][s][t]:
					people.append([h, t])
		
		# Para cada par de pessoas Hi e Hj 
		for i in people:
			Hi = i[0] # indice da pessoa i
			Xip = i[1] # indice do tempo de i no projeto p
			
			for j in people:
				Hj = j[0] # indice da pessoa j
				Xjp = j[1] # indice do tempo de j no projeto p
				
				#if Hi != Hj:
				e += S[Hi][Hj] * Xip * Xjp
		
		Rsp = 0
		
		# Para cada demanda de habilidade no projeto p
		for r in R[p]:
			Rsp += r
		
		# Calcula a eficiencia do projeto p e acumula na eficiencia global
		#efficiency += (((e / Rsp**2) + 1) / 2) * W[p]
		efficiency += (((e / Rsp**2) + 1) / 2) * (1/m)
		

	return efficiency
'''


# Retorna um array com a probabilidade de cada habilidade ser sorteada, baseada nas demandas pelas habilidades
def skill_probability(R, m, f):
	array_probability = []
	
	# Total de demandas de todas as habilidades
	total_req = 0
	
	for s in range(0, f):
		# Total de demandas da habilidade s
		s_req = 0
		
		for p in range(0, m):
			s_req += R[p][s]
		
		array_probability.append(s_req)
		
		total_req += s_req
	
	for i in range(0, len(array_probability)):
		array_probability[i] /= total_req
	
	return array_probability 


# Retorna um array com a probabilidade de cada habilidade ser sorteada, baseada nas demandas pelas habilidades
# VERSAO 2: coloca probabilidade 0 se uma habilidade possui demanda apenas em 1 projeto
def skill_probability_v2(R, m, f):
	array_probability = []
	
	# Total de demandas de todas as habilidades
	total_req = 0
	
	for s in range(0, f):
		# Total de demandas da habilidade s
		s_req = 0
		
		aux = 0
		
		for p in range(0, m):
			if R[p][s] != 0:
				s_req += R[p][s]
				
				aux += 1
		
		if aux == 1:
			s_req = 0
			
		array_probability.append(s_req)
		
		total_req += s_req
	
	for i in range(0, len(array_probability)):
		array_probability[i] /= total_req
	
	return array_probability 
	

# Retorna um array com a probabilidade de cada habilidade ser sorteada, baseada nas pessoas que possuem cada habilidade
def skill_people_probability(K, n, f):
	array_probability = []
	
	# Total de demandas de todas as habilidades
	total_req = 0
	
	for s in range(0, f):
		# Total de demandas da habilidade s
		s_req = 0
		
		for p in range(0, n):
			s_req += K[p][s]
		
		array_probability.append(s_req)
		
		total_req += s_req
	
	for i in range(0, len(array_probability)):
		array_probability[i] /= float(total_req)
	
	return array_probability 


# Dada uma habilidade skill, retorna um array com a probabilidade de cada projeto ser sorteado, baseada na demanda de skill
def project_probability(R, m, f, skill):
	array_probability = []
	
	# Total de demandas de todas as habilidades
	total_req = 0
	
	for p in range(0, m):
		# Demanda da habilidade skill no projeto p
		p_req = R[p][skill]
	
		array_probability.append(p_req)
	
		total_req += p_req
		
	for i in range(0, len(array_probability)):
		array_probability[i] /= total_req	
	
	return array_probability


# Dadas 2 habilidades, retorna um array com a probabilidade de cada projeto ser sorteado, baseada na soma das demandas das 
# 2 habilidades
def project_2_skill_probability(R, m, f, skill1, skill2):
	array_probability = []
	
	# Total de demandas de todas as habilidades
	total_req = 0
	
	for p in range(0, m):
		# Soma as demandas das habilidades no projeto p
		p_req = R[p][skill1] + R[p][skill2]
	
		array_probability.append(p_req)
	
		total_req += p_req
		
	for i in range(0, len(array_probability)):
		array_probability[i] /= total_req	
	
	return array_probability


# Retorna um array com as propocoes de tempos fracionados e de tempos inteiros das demandas em R
def fractional_time_probability(R):
	# Array de 2 posicoes, onde a primeira indica a proporcao de tempos fracionados e a segunda indica a proporcao de 
	# tempos inteiros
	array_probability = [0, 0]
	
	total = 0.0
	
	for p in range(0, len(R)):
		for s in range(0, len(R[p])):
			if R[p][s] != 0:
				# Se tempo fracionado
				if (R[p][s] % 1) != 0:
					array_probability[0] += 1
				else:
					array_probability[1] += 1
				
				total += 1
	
	for i in range(0, len(array_probability)):
		array_probability[i] /= total
				
	return array_probability


# Dados dois projetos e uma habilidade, retorna um array com a probabilidade de cada tempo ser sorteado, baseado na quantidade
# de pessoas que estao no subgrupo de cada tempo daquela habilidade em cada projeto
def time_probability(project1, project2, skill, solution):
	array_probability = []
	
	# Dicts de tempos dos projetos 1 e 2 para a habilidade skill
	time_proj1 = solution[project1][skill]
	time_proj2 = solution[project2][skill]
	
	# Lista de tempos
	times = list(solution[project1][skill])
	
	# Total de pessoas
	total_people = 0
	
	for t in times:
		# Quantidade de pessoas com o tempo t nos projetos 1 e 2 
		qt_people_proj1_t = len(time_proj1[t])
		qt_people_proj2_t = len(time_proj2[t])
		
		# Se nao houver pessoas no subgrupo do tempo t em um dos dois projetos, a probabilidade de escolher esse tempo sera 0
		if (qt_people_proj1_t == 0) or (qt_people_proj2_t == 0):
			array_probability.append(0)
			
		else:
			people_t = qt_people_proj1_t + qt_people_proj2_t

			array_probability.append(people_t)
			
			total_people += people_t
	
	if total_people != 0:
		for i in range(0, len(array_probability)):
			array_probability[i] /= total_people
		
	return array_probability


# Dados dois projetos e 2 habilidades, retorna um array com a probabilidade de cada tempo ser sorteado, baseado na quantidade
# de pessoas que estao no subgrupo de cada tempo das 2 habilidades em cada projeto
def time_2_skill_probability(project1, project2, skill1, skill2, solution):
	array_probability = []
	
	# Dicts de tempos dos projetos 1 e 2 para as 2 habilidades
	time_proj1_skill1 = solution[project1][skill1]
	time_proj1_skill2 = solution[project1][skill2]
	time_proj2_skill1 = solution[project2][skill1]
	time_proj2_skill2 = solution[project2][skill2]
	
	# Lista de tempos
	times = list(solution[project1][skill1])
	
	# Total de pessoas
	total_people = 0
	
	for t in times:
		# Quantidade de pessoas com o tempo t nos projetos 1 e 2 
		qt_people_proj1_t = len(time_proj1_skill1[t]) + len(time_proj1_skill2[t])
		qt_people_proj2_t = len(time_proj2_skill1[t]) + len(time_proj2_skill2[t])
		
		people_t = qt_people_proj1_t + qt_people_proj2_t

		array_probability.append(people_t)
		
		total_people += people_t
	
	if total_people != 0:
		for i in range(0, len(array_probability)):
			array_probability[i] /= total_people
		
	return array_probability


# Calcula a probabilidade de cada solucao ser selecionada dentro de uma populacao, de acordo com sua eficiencia
def solution_probability(population_quality):
	array_probability = []
	
	total = 0 
	
	for solution_quality in population_quality:
		array_probability.append(solution_quality)
		
		total += solution_quality
	
	for i in range(0, len(array_probability)):
		array_probability[i] /= total
	
	return array_probability


# Calcula a probabilidade de um projeto ser selecionado para sofrer mutacao; quanto menor a eficiência do projeto, mais alta
# a probabilidade.
def project_mutation_probability(solution, S, R, m):
	array_probability = []

	total = 0
	
	for project in solution:
		e = project_efficiency(solution, project, S, R)
		p = 1 - e

		array_probability.append(p)
		
		total += p
	
	if total != 0:
		for i in range(0, len(array_probability)):
			array_probability[i] /= total
	
	return array_probability


# Dado um array de probabilidades, sera feito um sorteio probabilistico e retornado um valor do array
def probabilistic_draw(array_probability):
	# Numero aleatorio entre 0.0 e 1.0
	prob = random.random()
	
	# Valor escolhido
	value = -1
	
	aux = 0
	
	# Valor escolhida de acordo com a probabilidade
	for v in range(0, len(array_probability)):
		prob_v = array_probability[v]
		
		if prob >= aux and prob < (aux + prob_v):
			# Escolhe o valor v
			value = v
			break
		
		aux += prob_v
	
	return value


# Verificar se o tempo de uma pessoa aumentou na solucao, para um dado projeto e habilidade
def verify_people_time(solution, project, skill, person):
	#print("###########################################")
	#print(solution)

	times = list(solution[project][skill])
	
	# Armazena os tempos em que a pessoa esta na solucao
	people_times = []
	
	for t in times:
		people = solution[project][skill][t]
		
		for p in people:
			if p == person:
				people_times.append(t)
	
	if len(people_times) > 1:
		# Probabilidade de somar os tempos
		prob = random.random()
		
		if prob < 0.5:
			#print("UNIU OS TEMPOS DA PESSOA " + str(person))
			time = sum(people_times)
			#print(people_times)
			
			# Aloco a pessoa no novo tempo
			solution[project][skill][time].append(person)
			
			# Removo a pessoa dos tempos antigos
			for t in people_times:
				solution[project][skill][t].remove(person)
	
	return solution


# Troca duas pessoas de equipes diferentes que estao exercendo a mesma habilidade na mesma fracao de tempo
def swap1(solution, R, m, f, D):
	#print("SOLUCAO RECEBIDA PELO SWAP 1")
	#print(solution)

	# Array de probabilidade das habilidades
	#skill_prob_array = skill_probability(R, m, f)
	skill_prob_array = skill_probability_v2(R, m, f)
	
	# Habilidade sorteada
	skill = probabilistic_draw(skill_prob_array)
	
	#print("Habilidade: " + str(skill))
	
	# Array de probabilidade dos projetos
	project_prob_array = project_probability(R, m, f, skill)
	
	#print("Prob projetos: " + str(project_prob_array))
	
	# Projeto 1 sorteado
	project1 = probabilistic_draw(project_prob_array)

	#print("Projeto: " + str(project1))	
	
	# Armazena a probabilidade do projeto 1
	prob_proj1 = project_prob_array[project1]
	
	# Altera a probabilidade do projeto 1 para 0, para nao ser sorteado de novo
	project_prob_array[project1] = 0
	
	cnt = 0
	
	# cnt armazena a quantidade de projetos com probabilidade diferente de 0
	for p in range(0, len(project_prob_array)):
		if project_prob_array[p] != 0.0:
			cnt += 1
	
	# Divide-se a probabilidade do projeto 1 entre os projetos que ja possuiam uma probabilidade maior que 0
	if cnt != 0:	
		prob = prob_proj1 / cnt
	
		for p in range(0, len(project_prob_array)):
			if project_prob_array[p] != 0.0:
				project_prob_array[p] += prob
	
		#print("Prob projetos: " + str(project_prob_array))	
		
		# Projeto 2 sorteado
		project2 = probabilistic_draw(project_prob_array)
		
		#print("Projeto: " + str(project2))
		
		# Array de probabilidade dos tempos
		time_prob_array = time_probability(project1, project2, skill, solution)

		#print("Prob tempos: " + str(time_prob_array))
		
		if sum(time_prob_array) != 0:
			# Tempo sorteado
			time_index = probabilistic_draw(time_prob_array)
			time = D[time_index]
		
			# Sorteia duas pessoas para a troca
			person1_index = random.randint(0, len(solution[project1][skill][time]) - 1)
			person2_index = random.randint(0, len(solution[project2][skill][time]) - 1)
	
			person1 = solution[project1][skill][time][person1_index]
			person2 = solution[project2][skill][time][person2_index]
			
			#print("Pessoas: " + str(person1) + " e " + str(person2))
			
			if person1 != person2:
				solution[project1][skill][time][person1_index] = person2 # Pessoa 2 vai para o projeto 1
				solution[project2][skill][time][person2_index] = person1 # Pessoa 1 vai para o projeto 2
				
				# Apos uma troca, o tempo de uma das duas pessoas pode aumentar, pois ela pode ter mais tempo alocado a uma equipe
				# em que ela ja estava trabalhando
				solution = verify_people_time(solution, project1, skill, person2)
				solution = verify_people_time(solution, project2, skill, person1)
				
				return solution
			
			'''	
			# Pessoas iguais sorteadas, o que fazer??
			else:
				print("SORTEIA OUTRA HABILIDADE?")
				print("VAI PRA OUTRO SWAP?")				
							
		# Nao ha 2 pessoas nos dois projetos com a habilidade escolhida e tempo igual
		else:
			print("SORTEIA OUTRA HABILIDADE?")
			print("VAI PRA OUTRO SWAP?")
	
	# Nao ha outro projeto com a habilidade sorteada
	else:
		print("SORTEIA OUTRA HABILIDADE?")
		print("VAI PRA OUTRO SWAP?")

		'''

# Troca duas pessoas que estao exercendo habilidades diferentes mas que possuem a habilidade que a outra pessoa esta exercendo 
# em equipes diferentes em frações de tempos iguais
def swap2(solution, R, m, f, D, K):
	# Array de probabilidade das habilidades
	skill_prob_array = skill_probability(R, m, f)
	
	# Habilidade 1 sorteada
	skill1 = probabilistic_draw(skill_prob_array)
	
	#print("Habilidade: " + str(skill1))
	
	# Armazena a probabilidade da habilidade 1
	prob_skill1 = skill_prob_array[skill1]
	
	# Altera a probabilidade da habilidade1 para 0, para nao ser sorteada de novo
	skill_prob_array[skill1] = 0
	
	cnt = 0
	
	# cnt armazena a quantidade de habilidades com probabilidade diferente de 0
	for s in range(0, len(skill_prob_array)):
		if skill_prob_array[s] != 0.0:
			cnt += 1
			
	# Distribui a probabilidade da habilidade1 igualmente entre as habilidades que ja possuiam probabilidade maior que 0
	if cnt != 0:		
		prob = prob_skill1 / cnt
		
		for s in range(0, len(skill_prob_array)):
			if skill_prob_array[s] != 0.0:
				skill_prob_array[s] += prob
		
		# Habilidade 2 sorteada
		skill2 = probabilistic_draw(skill_prob_array)
	
		#print("Habilidade: " + str(skill2))
		
		# Array de probabilidade dos projetos, baseado na demanda das 2 habilidades escolhidas
		project_prob_array = project_2_skill_probability(R, m, f, skill1, skill2)
		
		# Projeto 1 sorteado
		project1 = probabilistic_draw(project_prob_array)
		
		#print("Projeto: " + str(project1))
		
		# Armazena a probabilidade do projeto 1
		prob_proj1 = project_prob_array[project1]
		
		# Altera a probabilidade do projeto 1 para 0, para nao ser sorteado de novo
		project_prob_array[project1] = 0
		
		cnt = 0
		
		# cnt armazena a quantidade de projetos com probabilidade diferente de 0
		for p in range(0, len(project_prob_array)):
			if project_prob_array[p] != 0.0:
				cnt += 1
		
		# Divide-se a probabilidade do projeto 1 entre os projetos que ja possuiam uma probabilidade maior que 0
		if cnt != 0:	
			prob = prob_proj1 / cnt
		
			for p in range(0, len(project_prob_array)):
				if project_prob_array[p] != 0.0:
					project_prob_array[p] += prob
		
			# Projeto 2 sorteado
			project2 = probabilistic_draw(project_prob_array)
			
			#print("Projeto: " + str(project2))
			
			time_prob_array = time_2_skill_probability(project1, project2, skill1, skill2, solution)
			
			if sum(time_prob_array) != 0:
				# Tempo sorteado
				time_index = probabilistic_draw(time_prob_array)
				time = D[time_index]
				
				#print("Tempo: " + str(time))
				
				# pessoa_projeto_habilidadeExercida_possuiHabilidade
				people_project1_skill1 = solution[project1][skill1][time]; people_project1_skill1_skill2 = []
				people_project1_skill2 = solution[project1][skill2][time]; people_project1_skill2_skill1 = []
				people_project2_skill1 = solution[project2][skill1][time]; people_project2_skill1_skill2 = []
				people_project2_skill2 = solution[project2][skill2][time]; people_project2_skill2_skill1 = []
				
				# Armazena as pessoas do projeto 1, que estao exercendo a habilidade 1 e possuem a habilidade 2
				for p in people_project1_skill1:
					if K[p][skill2] != 0:
						people_project1_skill1_skill2.append(p)
				
				# Armazena as pessoas do projeto 1, que estao exercendo a habilidade 2 e possuem a habilidade 1
				for p in people_project1_skill2:
					if K[p][skill1] != 0:
						people_project1_skill2_skill1.append(p)
				
				# Armazena as pessoas do projeto 2, que estao exercendo a habilidade 1 e possuem a habilidade 2
				for p in people_project2_skill1:
					if K[p][skill2] != 0:
						people_project2_skill1_skill2.append(p)
				
				# Armazena as pessoas do projeto 2, que estao exercendo a habilidade 2 e possuem a habilidade 1
				for p in people_project2_skill2:
					if K[p][skill1] != 0:
						people_project2_skill2_skill1.append(p)
				
				# flag1 representa o primeiro sorteio possivel e flag2 representa o segundo sorteio possivel	
				flag1 = False; flag2 = False
				
				# Primeiro sorteio possivel: 
				# uma pessoa do projeto 1, exercendo a habilidade 1 e que possui a habilidade 2
				# uma pessoa do projeto 2, exercendo a habilidade 2 e que possui a habilidade 1
				if (len(people_project1_skill1_skill2) != 0) and (len(people_project2_skill2_skill1) != 0):
					flag1 = True
				
				# Segundo sorteio possivel: 
				# uma pessoa do projeto 1, exercendo a habilidade 2 e que possui a habilidade 1
				# uma pessoa do projeto 2, exercendo a habilidade 1 e que possui a habilidade 2
				elif (len(people_project1_skill2_skill1) != 0) and (len(people_project2_skill1_skill2) != 0):	
					flag2 = True
				
				
				person1 = -1
				person2 = -1
				
				if flag1 == True and flag2 == True:
					# 50% de probabilidade para cada sorteio
					prob = random.random()
					
					# Sorteio 1
					if prob < 0.5:
						person1_index = random.randint(0, len(people_project1_skill1_skill2) - 1)
						person1 = people_project1_skill1_skill2[person1_index]
						
						person2_index = random.randint(0, len(people_project2_skill2_skill1) - 1)
						person2 = people_project2_skill2_skill1[person2_index]
						
						#print("Pessoas: " + str(person1) + " e " + str(person2))
						
						if person1 != person2:
							solution[project1][skill1][time].append(person2) # Pessoa 2 vai para o projeto 1
							solution[project2][skill2][time].remove(person2) # Pessoa 2 é removida do projeto 2
							
							solution[project2][skill2][time].append(person1) # Pessoa 1 vai para o projeto 2
							solution[project1][skill1][time].remove(person1) # Pessoa 1 é removida do projeto 1
							
							#print("Pessoa " + str(person2) + "(Projeto " + str(project2) + " , Habilidade " + str(skill2) + ")" + " -> (Projeto " + str(project1) + ", Habilidade " + str(skill1) + ")")
							#print("Pessoa " + str(person1) + "(Projeto " + str(project1) + " , Habilidade " + str(skill1) + ")" + " -> (Projeto " + str(project2) + ", Habilidade " + str(skill2) + ")")
							
							# Apos uma troca, o tempo de uma das duas pessoas pode aumentar, pois ela pode ter mais tempo alocado 
							# a uma equipe em que ela ja estava trabalhando
							solution = verify_people_time(solution, project1, skill1, person2)
							solution = verify_people_time(solution, project2, skill2, person1)
							
							return solution
						
						# Pessoas sortedas sao iguais
						#else:
						#	print("O QUE FAZER?")
						
					# Sorteio 2
					else:
						person1_index = random.randint(0, len(people_project1_skill2_skill1) - 1)
						person1 = people_project1_skill2_skill1[person1_index]
						
						person2_index = random.randint(0, len(people_project2_skill1_skill2) - 1)
						person2 = people_project2_skill1_skill2[person2_index]
						
						#print("Pessoas: " + str(person1) + " e " + str(person2))
						
						if person1 != person2:						
							solution[project1][skill2][time].append(person2) # Pessoa 2 vai para o projeto 1
							solution[project2][skill1][time].remove(person2) # Pessoa 2 é removida do projeto 2
							
							solution[project2][skill1][time].append(person1) # Pessoa 1 vai para o projeto 2
							solution[project1][skill2][time].remove(person1) # Pessoa 1 é removida do projeto 1
							
							#print("Pessoa " + str(person2) + "(Projeto " + str(project2) + " , Habilidade " + str(skill1) + ")" + " -> (Projeto " + str(project1) + ", Habilidade " + str(skill2) + ")")
							#print("Pessoa " + str(person1) + "(Projeto " + str(project1) + " , Habilidade " + str(skill2) + ")" + " -> (Projeto " + str(project2) + ", Habilidade " + str(skill1) + ")")
							
							# Apos uma troca, o tempo de uma das duas pessoas pode aumentar, pois ela pode ter mais tempo alocado 
							# a uma equipe em que ela ja estava trabalhando
							solution = verify_people_time(solution, project1, skill2, person2)
							solution = verify_people_time(solution, project2, skill1, person1)
							
							return solution
						
						# Pessoas sortedas sao iguais
						#else:
						#	print("O QUE FAZER?")
							
				elif flag1 == True and flag2 == False:
					# Sorteio 1
					person1_index = random.randint(0, len(people_project1_skill1_skill2) - 1)
					person1 = people_project1_skill1_skill2[person1_index]
					
					person2_index = random.randint(0, len(people_project2_skill2_skill1) - 1)
					person2 = people_project2_skill2_skill1[person2_index]
					
					#print("Pessoas: " + str(person1) + " e " + str(person2))
					
					if person1 != person2:
						solution[project1][skill1][time].append(person2) # Pessoa 2 vai para o projeto 1
						solution[project2][skill2][time].remove(person2) # Pessoa 2 é removida do projeto 2
						
						solution[project2][skill2][time].append(person1) # Pessoa 1 vai para o projeto 2
						solution[project1][skill1][time].remove(person1) # Pessoa 1 é removida do projeto 1
							
						#print("Pessoa " + str(person2) + "(Projeto " + str(project2) + " , Habilidade " + str(skill2) + ")" + " -> (Projeto " + str(project1) + ", Habilidade " + str(skill1) + ")")
						#print("Pessoa " + str(person1) + "(Projeto " + str(project1) + " , Habilidade " + str(skill1) + ")" + " -> (Projeto " + str(project2) + ", Habilidade " + str(skill2) + ")")
						
						# Apos uma troca, o tempo de uma das duas pessoas pode aumentar, pois ela pode ter mais tempo alocado 
						# a uma equipe em que ela ja estava trabalhando
						solution = verify_people_time(solution, project1, skill1, person2)
						solution = verify_people_time(solution, project2, skill2, person1)
						
						return solution
					
					# Pessoas sortedas sao iguais
					#else:
						#print("O QUE FAZER?")
						
				elif flag1 == False and flag2 == True:
					# Sorteio 2
					person1_index = random.randint(0, len(people_project1_skill2_skill1) - 1)
					person1 = people_project1_skill2_skill1[person1_index]
					
					person2_index = random.randint(0, len(people_project2_skill1_skill2) - 1)
					person2 = people_project2_skill1_skill2[person2_index]
					
					#print("Pessoas: " + str(person1) + " e " + str(person2))
					
					if person1 != person2:
						solution[project1][skill2][time].append(person2) # Pessoa 2 vai para o projeto 1
						solution[project2][skill1][time].remove(person2) # Pessoa 2 é removida do projeto 2
						
						solution[project2][skill1][time].append(person1) # Pessoa 1 vai para o projeto 2
						solution[project1][skill2][time].remove(person1) # Pessoa 1 é removida do projeto 1
						
						#print("Pessoa " + str(person2) + "(Projeto " + str(project2) + " , Habilidade " + str(skill1) + ")" + " -> (Projeto " + str(project1) + ", Habilidade " + str(skill2) + ")")
						#print("Pessoa " + str(person1) + "(Projeto " + str(project1) + " , Habilidade " + str(skill2) + ")" + " -> (Projeto " + str(project2) + ", Habilidade " + str(skill1) + ")")
						
						# Apos uma troca, o tempo de uma das duas pessoas pode aumentar, pois ela pode ter mais tempo alocado 
						# a uma equipe em que ela ja estava trabalhando
						solution = verify_people_time(solution, project1, skill2, person2)
						solution = verify_people_time(solution, project2, skill1, person1)
						
						return solution
					'''
					# Pessoas sortedas sao iguais
					else:
						print("O QUE FAZER?")
						
				# Nao da para fazer a troca, pois nao existem tais pessoas 
				else:
					print("O QUE FAZER?")
						
			# Nao ha ninguem nos dois projetos com as duas habilidades (muito improvavel)	
			else:
				print("O QUE FAZER?")
			
		# Nao ha outro projeto com demandas das duas habilidades
		else:
			print("O QUE FAZER?")
	
	# Nao ha outra habilidade a ser sorteada (acho que nunca ocorre)
	else:
		print("O QUE FAZER?")
				'''


# Sao escolhidas 2 equipes e 2 pessoas que estao nessas 2 equipes com a mesma habilidade em cada equipe e os tempos delas 
# sao trocados, em alfa unidades
def swap3(solution, alfa, R, m, f):
	print('Swap 3 ocorreu')

	# Array de probabilidade das habilidades
	# ALTERAR PARA A VERSÃO 2 DO SKILL_PROBABILITY
	skill_prob_array = skill_probability(R, m, f)
	
	# Habilidade sorteada
	skill = probabilistic_draw(skill_prob_array)
	
	#print("Habilidade: " + str(skill))
	
	# Array de probabilidade dos projetos
	project_prob_array = project_probability(R, m, f, skill)
	
	#print("Prob projetos: " + str(project_prob_array))
	
	# Projeto 1 sorteado
	project1 = probabilistic_draw(project_prob_array)

	#print("Projeto: " + str(project1))	
	
	# Armazena a probabilidade do projeto 1
	prob_proj1 = project_prob_array[project1]
	
	# Altera a probabilidade do projeto 1 para 0, para nao ser sorteado de novo
	project_prob_array[project1] = 0
	
	cnt = 0
	
	# cnt armazena a quantidade de projetos com probabilidade diferente de 0
	for p in range(0, len(project_prob_array)):
		if project_prob_array[p] != 0.0:
			cnt += 1
	
	# Divide-se a probabilidade do projeto 1 entre os projetos que ja possuiam uma probabilidade maior que 0
	if cnt != 0:	
		prob = prob_proj1 / cnt
	
		for p in range(0, len(project_prob_array)):
			if project_prob_array[p] != 0.0:
				project_prob_array[p] += prob
	
		#print("Prob projetos: " + str(project_prob_array))	
		
		# Projeto 2 sorteado
		project2 = probabilistic_draw(project_prob_array)
		
		#print("Projeto: " + str(project2))
		
		# Armazena as pessoas do projeto 1 e do projeto 2, respectivamente, que estao exercendo a habilidade sorteada
		proj1_people = []
		proj2_people = []
		
		dict_people = dict()
		
		for t in solution[project1][skill]:
			for p in solution[project1][skill][t]:
				proj1_people.append(p)
				
				if (p in dict_people.keys()) == False:
					dict_people[p] = dict()
				
				if (project1 in dict_people[p].keys()) == False:
					dict_people[p][project1] = []
				
				dict_people[p][project1].append(t)
				
			for p in solution[project2][skill][t]:
				proj2_people.append(p)
				
				if (p in dict_people.keys()) == False:
					dict_people[p] = dict()
				
				if (project2 in dict_people[p].keys()) == False:
					dict_people[p][project2] = []
				
				dict_people[p][project2].append(t)
		
		proj1_people = set(proj1_people)
		
		# Seleciona somente as pessoas que estao em ambos os projetos, exercendo a mesma habilidade
		people = list(proj1_people.intersection(proj2_people))
		
		#print('Pessoas que posso trocar: ' + str(people))
		
		# Nao ha 2 pessoas possiveis para realizar a troca
		if len(people) < 2:
			#print("O QUE FAZER?")
			return
		
		person1_index = random.randint(0, len(people) - 1)
		person1 = people[person1_index]
		
		people.remove(person1)
		
		person2_index = random.randint(0, len(people) - 1)
		person2 = people[person2_index]
		
		#print("Pessoas: " + str(person1) + " e " + str(person2))
		#print(dict_people)
		time_person1_proj1_index = random.randint(0, len(dict_people[person1][project1]) - 1)
		time_person1_proj1 = dict_people[person1][project1][time_person1_proj1_index]
		
		time_person1_proj2_index = random.randint(0, len(dict_people[person1][project2]) - 1)
		time_person1_proj2 = dict_people[person1][project2][time_person1_proj2_index]
		
		time_person2_proj1_index = random.randint(0, len(dict_people[person2][project1]) - 1)
		time_person2_proj1 = dict_people[person2][project1][time_person2_proj1_index]
	
		time_person2_proj2_index = random.randint(0, len(dict_people[person2][project2]) - 1)
		time_person2_proj2 = dict_people[person2][project2][time_person2_proj2_index]
		
		new_time_person1_proj1 = time_person1_proj1 + alfa
		new_time_person2_proj1 = time_person2_proj1 - alfa
		
		new_time_person1_proj2 = time_person1_proj2 - alfa
		new_time_person2_proj2 = time_person2_proj2 + alfa
		
		solution[project1][skill][time_person1_proj1].remove(person1)
		solution[project1][skill][new_time_person1_proj1].append(person1)
		
		if new_time_person1_proj2 != 0:
			solution[project2][skill][time_person1_proj2].remove(person1)
			solution[project2][skill][new_time_person1_proj2].append(person1)
		else:
			solution[project2][skill][time_person1_proj2].remove(person1)
		
		if new_time_person2_proj1 != 0:
			solution[project1][skill][time_person2_proj1].remove(person2)
			solution[project1][skill][new_time_person2_proj1].append(person2)
		else:
			solution[project1][skill][time_person2_proj1].remove(person2)
		
		solution[project2][skill][time_person2_proj2].remove(person2)
		solution[project2][skill][new_time_person2_proj2].append(person2)		
		
		return solution
		
	# Nao ha outro projeto com demanda da habilidade sorteada
	#else:
		#print("O QUE FAZER?")


# Dadas 2 solucoes, escolhe-se um ponto de troca entre as equipes e gera duas novas solucoes a partir do crossover das 2 solucoes
# nos pontos de troca escolhidos. Verifica-se se ha alguma pessoa inviavel nesta nova solucao, ou seja, se esta trabalhando mais
# de 100% de tempo. Para cada pessoa inviavel, troca-se essa pessoa por uma pessoa de fora, que tenha tempo livre e a mesma
# habilidade (mutacao). 
def swap4(solution1, solution2, m, D, n, K):
	# Essa solucao pega a primeira parte da solucao 1 e a segunda parte da solucao 2
	new_solution1 = dict()
	
	# Essa solucao pega a primeira parte da solucao 2 e a segunda parte da solucao 1
	new_solution2 = dict()
	
	# Gera um ponto de crossover entre 0 e m - 2. Por exemplo, se temos 4 projetos (m = 4), entao temos (m - 1 = 3) pontos 
	# de troca. Mas como quero comecar do ponto 0, fazemos (m - 2 = 2), e teremos os possiveis pontos de troca 0, 1 e 2.
	point_crossover = random.randint(0, m - 2)
	
	# Quantidade de projetos que vem da parte 1
	proj_quant1 = 0
	
	# Quantidade de projetos que vem da parte 2
	proj_quant2 = 0
	
	aux1 = copy.deepcopy(solution1)
	aux2 = copy.deepcopy(solution2)
	
	# Pego a primeira parte (do comeco ate o ponto de crossover) da solucao 1 e da solucao 2
	for i in range(0, point_crossover + 1):
		#new_solution1[i] = solution1[i].copy()
		new_solution1[i] = aux1[i]
		
		#new_solution2[i] = solution2[i].copy()
		new_solution2[i] = aux2[i]
		
		proj_quant1 += 1
	
	# Pego a segunda parte (do ponto de crossover + 1 ate o fim) da solucao 2 e da solucao 1
	for i in range(point_crossover + 1, m):
		#new_solution1[i] = solution2[i].copy()
		new_solution1[i] = aux2[i]
		
		#new_solution2[i] = solution1[i].copy()
		new_solution2[i] = aux1[i]
		
		proj_quant2 += 1
	
	#print("####################### Solução 1 #############################")
	#print(solution1)

	#print("####################### Solução 2 #############################")
	#print(solution2)
	
	#print('####################### PONTO DE TROCA = ' + str(point_crossover) + ' #############################')

	#print("####################### Nova Solução 1 ##########################")
	#print(new_solution1)
	
	#print("####################### Nova Solução 2 ##########################")
	#print(new_solution2)
		
	people_time = dict()
	
	# Analiso a viabilidade dos individuos da parte 1 (por ser menor)
	if proj_quant1 <= proj_quant2:
		
		#print("-------------- Parte 1 <= Parte 2 --------------")
		
		#################### ANALISANDO A VIABILIDADE DA NOVA SOLUCAO 1
		#print("####################### Analisando a viabilidade da nova solução 1 ##########################")
		
		# Armazeno os tempos das pessoas na parte 1 da solucao
		for project_before_point in range(0, point_crossover + 1):
			for skill in new_solution1[project_before_point]:
				for time in new_solution1[project_before_point][skill]:
					for person in new_solution1[project_before_point][skill][time]:
						if (person in people_time.keys()) == False:
							people_time[person] = 0
						
						people_time[person] += time
		
		# Verifico se alguma pessoa da parte 1 esta tambem na parte 2, e se estiver, verifico se o tempo ultrapassa 100%
		for project_after_point in range(point_crossover + 1, m):
			for skill in new_solution1[project_after_point]:
				for time in new_solution1[project_after_point][skill]:
					index = 0
					
					for person in new_solution1[project_after_point][skill][time]:
						if (person in people_time.keys()):
							people_time[person] += time
							
							# Se o tempo da pessoa ultrapassar 100%, esta pessoa esta inviavel, e ocorrera uma mutacao nela
							if people_time[person] > 1:
								#print("Pessoa " + str(person) + " inviavel. Tempo: " + str(people_time[person]))
							
								new_solution1 = mutation(new_solution1, project_after_point, skill, time, person, D, n, K, index)
								
								# A mutacao deu certo							
								if new_solution1 != None:
									# Como a pessoa foi trocada, decrementa o tempo dela
									people_time[person] -= time
									
									#print("################### Nova Solução 1 - Depois da Mutação ###############")
									#print(new_solution1)
								
								# A mutacao nao deu certo, e a solucao continua inviavel
								else: 
									# Simplemente retorna None na solucao, ou permite a inviabialidade?
									#print("O QUE FAZER?")
									break
						
						#print("################### Nova Solução 1 - O Q CONTESENO ###############")
						#print(new_solution1)
						
						index += 1

					if new_solution1 == None:
						break
				
				if new_solution1 == None:
						break
			
			if new_solution1 == None:
						break

		
		#################### ANALISANDO A VIABILIDADE DA NOVA SOLUCAO 2
		#print("####################### Analisando a viabilidade da nova solução 2 ##########################")
		
		people_time = dict()
		
		# Armazeno os tempos das pessoas na parte 1 da solucao
		for project_before_point in range(0, point_crossover + 1):
			for skill in new_solution2[project_before_point]:
				for time in new_solution2[project_before_point][skill]:
					for person in new_solution2[project_before_point][skill][time]:
						if (person in people_time.keys()) == False:
							people_time[person] = 0
						
						people_time[person] += time
		
		# Verifico se alguma pessoa da parte 1 esta tambem na parte 2, e se estiver, verifico se o tempo ultrapassa 100%
		for project_after_point in range(point_crossover + 1, m):
			for skill in new_solution2[project_after_point]:
				for time in new_solution2[project_after_point][skill]:
					index = 0
				
					for person in new_solution2[project_after_point][skill][time]:
						if (person in people_time.keys()):
							people_time[person] += time
							
							# Se o tempo da pessoa ultrapassar 100%, esta pessoa esta inviavel, e ocorrera uma mutacao nela
							if people_time[person] > 1:
								#print("Pessoa " + str(person) + " inviavel. Tempo: " + str(people_time[person]))
							
								new_solution2 = mutation(new_solution2, project_after_point, skill, time, person, D, n, K, index)
								
								# A mutacao deu certo							
								if new_solution2 != None:
									# Como a pessoa foi trocada, decrementa o tempo dela
									people_time[person] -= time
									
									#print("################### Nova Solução 2 - Depois da Mutação ###############")
									#print(new_solution2)
								
								# A mutacao nao deu certo, e a solucao continua inviavel
								else: 
									# Simplemente retorna None na solucao, ou permite a inviabialidade?
									#print("O QUE FAZER?")
									break
									
						#print("################### Nova Solução 2 - O Q CONTESENO ###############")
						#print(new_solution2)
						
						index += 1

					if new_solution2 == None:
						break
				
				if new_solution2 == None:
						break
			
			if new_solution2 == None:
						break

	# Analiso a viabilidade dos individuos da parte 2 (por ser menor)
	else:
		#print("-------------- Parte 2 < Parte 1 --------------")
		
		#################### ANALISANDO A VIABILIDADE DA NOVA SOLUCAO 1
		#print("####################### Analisando a viabilidade da nova solução 1 ##########################")
		
		people_time = dict()
		
		# Armazeno os tempos das pessoas na parte 2 da solucao
		for project_after_point in range(point_crossover + 1, m):
			for skill in new_solution1[project_after_point]:
				for time in new_solution1[project_after_point][skill]:
					for person in new_solution1[project_after_point][skill][time]:
						if (person in people_time.keys()) == False:
							people_time[person] = 0
						
						people_time[person] += time
		
		# Verifico se alguma pessoa da parte 2 esta tambem na parte 1, e se estiver, verifico se o tempo ultrapassa 100%
		for project_before_point in range(0, point_crossover + 1):
			for skill in new_solution1[project_before_point]:
				for time in new_solution1[project_before_point][skill]:
					index = 0
				
					for person in new_solution1[project_before_point][skill][time]:
						if (person in people_time.keys()):
							people_time[person] += time
							
							# Se o tempo da pessoa ultrapassar 100%, esta pessoa esta inviavel, e ocorrera uma mutacao nela
							if people_time[person] > 1:
								#print("Pessoa " + str(person) + " inviavel. Tempo: " + str(people_time[person]))

								new_solution1 = mutation(new_solution1, project_before_point, skill, time, person, D, n, K, index)
								
								# A mutacao deu certo							
								if new_solution1 != None:
									# Como a pessoa foi trocada, decrementa o tempo dela
									people_time[person] -= time
									
									#print("################### Nova Solução 1 - Depois da Mutação ###############")
									#print(new_solution1)
								
								# A mutacao nao deu certo, e a solucao continua inviavel
								else: 
									# Simplemente retorna None na solucao, ou permite a inviabialidade?
									#print("O QUE FAZER?")
									break
						
						#print("################### Nova Solução 1 - O Q CONTESENO ###############")
						#print(new_solution1)
						
						index += 1

					if new_solution1 == None:
						break
				
				if new_solution1 == None:
						break
			
			if new_solution1 == None:
						break

		#################### ANALISANDO A VIABILIDADE DA NOVA SOLUCAO 2
		#print("####################### Analisando a viabilidade da nova solução 2 ##########################")
		
		people_time = dict()
		mutation_people = []
		
		# Armazeno os tempos as pessoas na parte 2 da solucao
		for project_after_point in range(point_crossover + 1, m):
			for skill in new_solution2[project_after_point]:
				for time in new_solution2[project_after_point][skill]:
					for person in new_solution2[project_after_point][skill][time]:
						if (person in people_time.keys()) == False:
							people_time[person] = 0
						
						people_time[person] += time
		
		# Verifico se alguma pessoa da parte 2 esta tambem na parte 1, e se estiver, verifico se o tempo ultrapassa 100%
		for project_before_point in range(0, point_crossover + 1):
			for skill in new_solution2[project_before_point]:
				for time in new_solution2[project_before_point][skill]:
					index = 0
				
					for person in new_solution2[project_before_point][skill][time]:
						if (person in people_time.keys()):
							people_time[person] += time
							
							# Se o tempo da pessoa ultrapassar 100%, esta pessoa esta inviavel, e ocorrera uma mutacao nela
							if people_time[person] > 1:
								#print("Pessoa " + str(person) + " inviavel. Tempo: " + str(people_time[person]))
							
								new_solution2 = mutation(new_solution2, project_before_point, skill, time, person, D, n, K, index)
								
								# A mutacao deu certo							
								if new_solution2 != None:
									# Como a pessoa foi trocada, decrementa o tempo dela
									people_time[person] -= time
									
									#print("################### Nova Solução 2 - Depois da Mutação ###############")
									#print(new_solution2)
								
								# A mutacao nao deu certo, e a solucao continua inviavel
								else: 
									# Simplemente retorna None na solucao, ou permite a inviabialidade?
									#print("O QUE FAZER?")
									break
									
						#print("################### Nova Solução 2 - O Q CONTESENO ###############")
						#print(new_solution2)
						
						index += 1

					if new_solution2 == None:
						break
				
				if new_solution2 == None:
						break
			
			if new_solution2 == None:
						break

	
	return new_solution1, new_solution2


# Dada uma solucao, proj, habilidade, tempo e uma pessoa que esta inviavel (com tempo maior que 1), realiza uma mutacao, ou seja, 
# troca essa pessoa inviavel por outra pessoa que tenha tempo disponivel e que possua esta habilidade
def mutation(solution, project, skill, time, person, D, n, K, index):
	dict_people_free_time = unassigned_people(solution, D, n)
	
	#print("Tempo Livre: " + str(dict_people_free_time))
	
	person2 = -1
	
	possible_people = []
	
	for t in dict_people_free_time:
		# Verifico apenas os tempos maiores ou iguais ao tempo da pessoa inviavel
		if t >= time:
		
			for p in dict_people_free_time[t]:
				# Se a pessoa p possui a mesma habilidade da pessoa inviavel, entao ela sera uma possivel pessoa para a troca
				if K[p][skill] != 0:
					person2 = p
					break
					#possible_people.append(p)
			
			if person2 != -1:
				break
	
	if person2 != -1:
	#if len(possible_people) != 0:
		#person2_index = random.randint(0, len(possible_people) - 1)
		#person2 = possible_people[person2_index]
	
		# Adiciona a nova pessoa na posicao da pessoa inviavel
		solution[project][skill][time][index] = person2
		
		#print("Pessoa " + str(person) + " removida.")
		#print("Pessoa " + str(person2) + " adicionada.")		
			
		return solution
			
	else:
		# Caso nao exista uma possivel pessoa para realizar a troca, retornar None
		return None


# Retorna o total de pessoas demandadas em todos os projetos, pegando sempre o teto (ceil)
def total_requirements(R):
	total = 0
	
	for project in range(len(R)):
		for skill in range(len(R[project])):
			total += math.ceil(R[project][skill])

	return total

# Define a probabilidade de cada swap ser sorteado, de acordo com a configuracao da instancia 
def swaps_probability(R, K, D, m, f, n):
	# Se a instancia nao é de multiplas habilidades (aux=1), se é (aux=0)
	aux = 1
	
	# Se a instancia possui tempo fracionado
	aux_fractional_time = 0
	
	if len(D) >= 2:
		aux_fractional_time = 1
		

	# Array de 4 posicoes, para os 4 swaps, respectivamente e a mutacao
	# Swap 1 comeca com 100%
	array_probability = [1, 0, 0, 0]

	################## Escolhendo a probabilidade do swap 2 ocorrer (posicao 1) #############################################
	skill_people_prob_array = skill_people_probability(K, n, f)
	
	w = 0
	s = 0 # sentinela
	
	if aux == 1:
		w = 0
	
	else:
		for skill_prob in skill_people_prob_array:
			# Se houver uma probabilidade maior ou igual a 70%, as habilidades nao estao muito distribuidas
			if skill_prob >= 0.7:
				w = 0.10
				s = 1
				
				break
		
		# Se s == 0, quer dizer que nao ha uma probabilidade maior ou igual a 70%, entao as habilidades estao bem distribuidas
		if s == 0:
			w = 0.20
			
			# Seleciono as duas habilidades de maior probabilidade, de acordo com as pessoas que as possuem
			copy_skill_array = skill_people_prob_array[:]
			
			skill1_prob = max(copy_skill_array)
			copy_skill_array.remove(skill1_prob)
			skill2_prob = max(copy_skill_array)
			
			skill1 = -1
			skill2 = -1
			
			for i in range(0, len(skill_people_prob_array)):
				if skill1 == -1 and skill_people_prob_array[i] == skill1_prob:
					skill1 = i
				
				elif skill2 == -1 and skill_people_prob_array[i] == skill2_prob:
					skill2 = i	
			
			project_2_skill_prob_array = project_2_skill_probability(R, m, f, skill1, skill2)
			
			for proj_prob in project_2_skill_prob_array:
				# Se houver uma probabilidade maior ou igual a 70%, as demandas das habilidades nao estao muito distribuidas
				if proj_prob >= 0.7:
					w += 0.10
			
			# Nao ha uma probabilidade maior ou igual a 70%, entao as demandas das habilidades estao bem distribuidas
			if w == 0.20:
				w += 0.20
		
		# Decremento w unidades da posicao 0, e incremento na posicao 1
		array_probability[0] -= w		
		array_probability[1] += w
		
	#########################################################################################################################
	
	################## Escolhendo a probabilidade do swap 3 ocorrer (posicao 2) #############################################
	'''
	w = 0
	
	fractional_time_array = fractional_time_probability(R)
	fractional_time = fractional_time_array[0]
	
	if aux_fractional_time == 0 or len(D) >= 2:
		array_probability[2] = 0
	
	else:
		# Se tem muito tempo fracionado
		if fractional_time >= 0.40:
			w = 0.05
		else:
			w = 0.015
		
		# Se cada pessoa possui apenas 1 habilidade, esse swap fica melhor
		if aux == 1:
			w += 0.035		
		
		if aux != 1:
			# Decremento w unidades das posicoes 0 e 1, e incremento 2w na posicao 2
			array_probability[0] -= w
			array_probability[1] -= w
		
		else:
			array_probability[0] -= (2*w)
		
		array_probability[2] += (2*w)
	'''
	#########################################################################################################################
	
	
	################## Escolhendo a probabilidade do swap 4 ocorrer (posicao 3) #############################################
	'''
	w = 0
	
	# Se tem pouco tempo fracionado
	if fractional_time < 0.40:
		w = 0.05
	else:
		w = 0.015		
	
	if aux != 1:
		# O swap 1 fica muito ruim quando temos D = [0.25, 0.5, 0.75, 1.0], entao vamos aumentar a prob do swap 4 e 
		# diminuir do swap 1
		if s == 0:
			''' '''
			if len(D) == 4:
				w += 0.10
				
			# Se nao possui tempo fracionado, ou seja, D = [0,1], pode ser melhor aumentar a probabilidade do swap 4
			if aux_fractional_time == 0:
				w += 0.10
			''' '''
			
			w += 0.15
			
		# Decremento w unidades das posicoes 0 e 1, e incremento na posicao 3
		array_probability[0] -= w
		array_probability[1] -= w
	
	else:
		
		# Se nao possui tempo fracionado, ou seja, D = [0,1], pode ser melhor aumentar a probabilidade do swap 4
		if aux_fractional_time == 0:
			w += 0.10	
		
		if len(D) == 1:
			w += 0.20
		
		# O swap 1 fica muito ruim quando temos D = [0.25, 0.5, 0.75, 1.0], entao vamos aumentar a prob do swap 4 e 
		# diminuir do swap 1
		if len(D) == 4:
			w += 0.485
			#w += 0.385
			#w += 0.30
		
		if len(D) == 2:
			w += 0.385
			
		array_probability[0] -= (2*w)
	
	
	array_probability[3] += (2*w)
	'''
	#########################################################################################################################
	
	################## Escolhendo a probabilidade da mutacao (posicao 4) #############################################
	
	
	return array_probability


# Retorna um dict com o tempo livre de cada pessoa, de acordo com uma solucao
def unassigned_people(solution, D, n):
	# Este dict armazena quanto tempo livre tem cada pessoa; sendo agrupado por tempo, assim: {0: [Pessoas], 0.25: [], ..., 1: []}
	dict_people_free_time = dict()
	
	# Dict auxiliar do tipo {Pessoa: Tempo}
	dict_aux = dict()
	
	dict_people_free_time[0] = []
	
	for time in D:
		dict_people_free_time[time] = []
	
	# Todas as pessoas primeiramente estao com tempo 100%
	for person in range(0, n):
		dict_people_free_time[1].append(person)
		
		dict_aux[person] = 1
	
	# De acordo com o tempo alocado na solucao, decrementa do tempo da pessoa
	for project in solution:
		for skill in solution[project]:
			for time in solution[project][skill]:
				for person in solution[project][skill][time]:
					# Recupero o tempo restante (antigo) da pessoa
					time_person = dict_aux[person]
					
					# Tempo restante (novo) da pessoa 
					time_left = time_person - time
					
					# No swap4, ao unirmos 2 solucoes, alguma pessoa pode ficar inviavel. Nesse caso, o time_left dela ficaria
					# negativo. Para fins de organizacao, setamos seu time_left para 0, uma vez que ela realmente nao tera
					# tempo livre. Ocorrera mutacao nessa pessoa. 
					if time_left < 0:
						time_left = 0
					
					# Remove a pessoa do tempo antigo
					dict_people_free_time[time_person].remove(person)					
					
					# Adiciona a pessoa ao novo tempo restante
					dict_people_free_time[time_left].append(person)
					
					# Atualiza o tempo restante da pessoa no dict_aux
					dict_aux[person] = time_left
	
	return dict_people_free_time


# Verifica se uma solucao eh viavel
def feasible(solution, R):
	#### Verificando se cada pessoa nao esta alocada em mais de 100% de tempo ###
	
	people_time = dict()

	for project in solution:
		for skill in solution[project]:
			demand = 0

			for time in solution[project][skill]:
				for person in solution[project][skill][time]:
					if (person in people_time.keys()) == False:
						people_time[person] = 0
					
					people_time[person] += time
					demand += time
					
			if demand != R[project][skill]:
				print ("Inviável por demanda!")
				return False
	
	for person in people_time:
		time = people_time[person]
		
		if time > 1:
			#print("------------------------------------------------------------------------------------------------------")
			#print(solution)
			print("Inviável por tempo!")
			return False
			
	return True


# Le as instancias
def read_data(S, K, R, D):
	'''
	tupS = read_S('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instanciaTeste/6Vertices/A1/S.txt')
	tupK = read_K('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instanciaTeste/6Vertices/A1/K.txt')
	tupR = read_R('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instanciaTeste/6Vertices/A1/R.txt')
	tupD = read_D('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instanciaTeste/6Vertices/A1/D.txt')
	tupW = read_W('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instanciaTeste/6Vertices/A1/W.txt')
	'''
	
	'''
	tupS = read_S('/home/marilia/Dropbox/TCC - AGRUPAMENTO/CLAIO - Modelo do artigo/ModelagemLinear/instancias/50Vertices/Classe5/D/II/S.txt')
	tupK = read_K('/home/marilia/Dropbox/TCC - AGRUPAMENTO/CLAIO - Modelo do artigo/ModelagemLinear/instancias/50Vertices/Classe5/D/II/K.txt')
	tupR = read_R('/home/marilia/Dropbox/TCC - AGRUPAMENTO/CLAIO - Modelo do artigo/ModelagemLinear/instancias/50Vertices/Classe5/D/II/R.txt')
	tupD = read_D('/home/marilia/Dropbox/TCC - AGRUPAMENTO/CLAIO - Modelo do artigo/ModelagemLinear/instancias/50Vertices/Classe5/D/II/D.txt')
	
	optimal = 0.988889
	'''
	
	optimal = 0
	
	tupS = read_S(S)
	tupK = read_K(K)
	tupR = read_R(R)
	tupD = read_D(D)
	
	S, n = tupS
	K, f = tupK
	R, m = tupR
	D, t, alfa = tupD
	#W, p = tupW
	
	#return (S, K, R, D, n, f, m, t, alfa, W, p)
	return (S, K, R, D, n, f, m, t, alfa), optimal
	
	
def genetic_algorithm(data, optimal):
	#S, K, R, D, n, f, m, t, alfa, W, p = data
	S, K, R, D, n, f, m, t, alfa = data
	
	# Gera as probabilidades de cada swap, de acordo com a instancia
	array_swaps_prob = swaps_probability(R, K, D, m, f, n)
	
	# Gera a probabilidade de mutacao, baseada no numero de pessoas ociosas
	prob_mut = 0
	
	free_people = n - total_requirements(R)
	
	if free_people == 0:
		prob_mutation = 0
		
	elif free_people > 10:
		prob_mutation = 0.20
		
	else:
		prob_mutation = 0.10
	
	if len(D) == 1:
		prob_mutation += 0.10
	
	#print(array_swaps_prob)
	#print(prob_mutation)
	#return None, None
	
	swap1_err = 0; qt_swap1 = 0
	swap2_err = 0; qt_swap2 = 0
	swap3_err = 0; qt_swap3 = 0
	swap4_err = 0; qt_swap4 = 0
	
	# Gera as solucoes iniciais. No momento, sao 50
	solutions = initials_solutions(S, K, R, D, n, f, m, t, alfa)
		
	new_population = solutions
	
	#x = 0
	
	best_solution_quality = 0
	conv = 0
	
	for i in range(0, 1000):
	#while x < 1000:
		population = new_population
		new_population = []
		
		if conv >= 200:
			break
		
		population_quality = []
		
		# Calcula a aptidao das solucoes da populacao
		for solution in population:
			#quality = global_efficiency(solution, S, R, alfa, W)
			quality = global_efficiency(solution, S, R, m)
			
			population_quality.append(quality)
			
			'''
			# Se solucao otima encontrada
			if quality == optimal:
				print("Encontrou o ótimo! :D")
				return solution, quality, swap1_err, swap2_err, swap3_err, swap4_err, array_swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4
			'''
		
		# Gero o array de probabilidades de cada solucao
		solution_prob_array = solution_probability(population_quality)
		
		# Sorteio um swap; O sorteio gera a partir do valor 0, que significa o swap 1, e assim por diante. 
		# Por isso o +1 no final.
		swap = probabilistic_draw(array_swaps_prob) + 1
		
		new_sol = None
		new_sol1 = None
		new_sol2 = None
		
		if swap == 1:
			qt_swap1 += 1
		
			#solution_index = random.randint(0, len(population) - 1)
			solution_index = probabilistic_draw(solution_prob_array)
			solution = population[solution_index]
			
			new_sol = swap1(solution, R, m, f, D)
			
			if new_sol == None:
				swap1_err += 1
			
			if ((new_sol != None) and ((feasible(new_sol, R)) == False)):
				print("---------------------------------SOLUCAO INVIAVEL GERADA PELO SWAP 1")
				print(new_sol)
				return
		
		elif swap == 2:
			qt_swap2 += 1
		
			#solution_index = random.randint(0, len(population) - 1)
			solution_index = probabilistic_draw(solution_prob_array)
			solution = population[solution_index]
		
			new_sol = swap2(solution, R, m, f, D, K)
			
			if new_sol == None:
				swap2_err += 1
			
			if ((new_sol != None) and ((feasible(new_sol, R)) == False)):
				print("---------------------------------SOLUCAO INVIAVEL GERADA PELO SWAP 2")
				print(new_sol)
				return
		
		elif swap == 3:
			qt_swap3 += 1
			
			#solution_index = random.randint(0, len(population) - 1)
			solution_index = probabilistic_draw(solution_prob_array)
			solution = population[solution_index]
			
			new_sol = swap3(solution, alfa, R, m, f)
			
			if new_sol == None:
				swap3_err += 1
			
			if ((new_sol != None) and ((feasible(new_sol, R)) == False)):
				print("---------------------------------SOLUCAO INVIAVEL GERADA PELO SWAP 3")
				print(new_sol)
				return
		
		else:
			qt_swap4 += 1
		
			#solution_index = random.randint(0, len(population) - 1)
			solution_index = probabilistic_draw(solution_prob_array)
			solution = population[solution_index]
			
			# Armazena a probabilidade da solucao 1
			prob_solution = solution_prob_array[solution_index]
			
			# Altera a probabilidade da solucao 1 para 0, para nao ser sorteado de novo
			solution_prob_array[solution_index] = 0
			
			cnt = 0
			
			# cnt armazena a quantidade de solucoes com probabilidade diferente de 0
			for p in range(0, len(solution_prob_array)):
				if solution_prob_array[p] != 0.0:
					cnt += 1
			
			# Divide-se a probabilidade da solucao 1 entre as solucoes que ja possuiam uma probabilidade maior que 0
			if cnt != 0:	
				prob = prob_solution / cnt
			
				for p in range(0, len(solution_prob_array)):
					if solution_prob_array[p] != 0.0:
						solution_prob_array[p] += prob
			
			#solution2_index = random.randint(0, len(population) - 1)
			solution2_index = probabilistic_draw(solution_prob_array)
			solution2 = population[solution2_index]
			
			#print("-------------- SOLUCAO 1 ---------------")
			#print(solution)
			
			#print("-------------- SOLUCAO 2 ---------------")
			#print(solution2)
			
			new_sol1, new_sol2 = swap4(solution, solution2, m, D, n, K)
			
			if new_sol1 == None:
				swap4_err += 1
			
			if new_sol2 == None:
				swap4_err += 1
			
			if ((new_sol1 != None) and ((feasible(new_sol1, R)) == False)):
				print("---------------------------------SOLUCAO 1 INVIAVEL GERADA PELO SWAP 4")
				print(new_sol1)
				return
			
			if ((new_sol2 != None) and ((feasible(new_sol2, R)) == False)):
				print("---------------------------------SOLUCAO 2 INVIAVEL GERADA PELO SWAP 4")
				print(new_sol2)
				return
			
			if feasible(solution) == False:
				print("ACHEEEEEEEEEEEEEEEEEEEEEEEI 1")
				print(solution)
				return
			
			if feasible(solution2) == False:
				print("ACHEEEEEEEEEEEEEEEEEEEEEEEI 2")
				print(solution2)
				return
				
		
		'''		
		if swap != 4 and new_sol == None:
			x -= 1
		
		if swap == 4 and new_sol1 == None and new_sol2 == None:
			x -= 1
		'''
		
		if swap != 4 and new_sol != None:
			#x += 1
			#print(new_sol)
			#q = global_efficiency(new_sol, S, R, alfa, W)
			#q = global_efficiency(new_sol, S, R)
			#print("Qualidade: " + str(q))		
			
			prob = random.random()
			
			if prob <= prob_mutation:
				proj_prob = project_mutation_probability(new_sol, S, R, m)
				
				if sum(proj_prob) != 0:
					project = probabilistic_draw(proj_prob)
					#project = random.randint(0, m-1)
					
					skill = random.randint(0, f-1)
					time_index = random.randint(0, t-1)
					time = D[time_index]
					
					if len(new_sol[project][skill][time]) > 0:
						person_index = random.randint(0, len(new_sol[project][skill][time]) - 1)
						person = new_sol[project][skill][time][person_index]
					
						aux = mutation(new_sol, project, skill, time, person, D, n, K, person_index)
						
						if aux != None:
							new_sol = aux
							
							if ((feasible(new_sol, R)) == False):
								print("---------------------------------SOLUCAO INVIAVEL GERADA PELA MUTACAO")
								print(new_sol)
								return
				
			population.append(new_sol)
		
		if swap == 4:
			if new_sol1 != None:
				#x += 1
			
				#print("SOLUCAO 1 GERADA PELO SWAP 4")
				#print(new_sol1)
				#q = global_efficiency(new_sol1, S, R, alfa, W)
				#q = global_efficiency(new_sol1, S, R)
				#print("Qualidade: " + str(q))
				
				prob = random.random()
			
				if prob <= prob_mutation:
					proj_prob = project_mutation_probability(new_sol1, S, R, m)
					
					if sum(proj_prob) != 0:
						project = probabilistic_draw(proj_prob)
						#project = random.randint(0, m-1)
						
						skill = random.randint(0, f-1)
						time_index = random.randint(0, t-1)
						time = D[time_index]
						
						if len(new_sol1[project][skill][time]) > 0:
							person_index = random.randint(0, len(new_sol1[project][skill][time]) - 1)
							person = new_sol1[project][skill][time][person_index]
						
							aux = mutation(new_sol1, project, skill, time, person, D, n, K, person_index)
							
							if aux != None:
								new_sol1 = aux
								
								if ((feasible(new_sol1, R)) == False):
									print("---------------------------------SOLUCAO INVIAVEL GERADA PELA MUTACAO")
									print(new_sol1)
									return
				
				population.append(new_sol1)
			
			if new_sol2 != None:
				#x += 1
			
				#print("SOLUCAO 2 GERADA PELO SWAP 4")
				#print(new_sol2)
				#q = global_efficiency(new_sol2, S, R, alfa, W)
				#q = global_efficiency(new_sol2, S, R)
				#print("Qualidade: " + str(q))
				
				prob = random.random()
				
				if prob <= prob_mutation:
					proj_prob = project_mutation_probability(new_sol2, S, R, m)
					
					if sum(proj_prob) != 0:
						project = probabilistic_draw(proj_prob)
						#project = random.randint(0, m-1)
						
						skill = random.randint(0, f-1)
						time_index = random.randint(0, t-1)
						time = D[time_index]
						
						if len(new_sol2[project][skill][time]) > 0:
							person_index = random.randint(0, len(new_sol2[project][skill][time]) - 1)
							person = new_sol2[project][skill][time][person_index]
						
							aux = mutation(new_sol2, project, skill, time, person, D, n, K, person_index)
							
							if aux != None:
								new_sol2 = aux
								
								if ((feasible(new_sol2, R)) == False):
									print("---------------------------------SOLUCAO INVIAVEL GERADA PELA MUTACAO")
									print(new_sol2)
									return
				
				population.append(new_sol2)
		
	
		solution_quality = [x for x in range(len(population))]
		i = 0
		
		for solution in population:
			#q = global_efficiency(solution, S, R, alfa, W)
			q = global_efficiency(solution, S, R, m)
		
			solution_quality[i] = q	
			
			i += 1
			
		#j = 0
		
		aux_list = copy.copy(solution_quality)
		aux_list = sorted(aux_list, reverse = True)
		
		# ESCOLHER APENAS OS X (45) MELHORES (E ALGUNS RUINS?) PARA A PROXIMA POPULACAO
		best = aux_list[:45]
		
		if len(aux_list) - 10 != 0:
			worse = aux_list[-5:]
		
		# Se o melhor individuo dessa populacao for igual ao melhor individuo já encontrado, aumento conv (pois pode chegar em convergncia)
		if best_solution_quality == best[0]:
			conv += 1
		
		# Se o melhor individuo dessa populacao for maior que o ultimo melhor individuo, entao este passa a ser o novo melhor e conv = 0, para recomecar a contagem da convergencia
		elif best_solution_quality < best[0]:
			best_solution_quality = best[0]
			conv = 0
		
		for i in best:
			solution_index = solution_quality.index(i)
			solution = population[solution_index]
			
			new_population.append(solution)
			
			#j += 1	
			
			#if j == 10:
			#	break
		
		for i in worse:
			solution_index = solution_quality.index(i)
			solution = population[solution_index]
			
			new_population.append(solution)
		

		
	# Retorna a melhor solucao		
	max_quality = 0
	max_solution = None
	
	for solution in population:
		#q = global_efficiency(solution, S, R, alfa, W)
		q = global_efficiency(solution, S, R, m)
	
		if q > max_quality:
			max_quality = q
			max_solution = solution
	
	
	return max_solution, max_quality, swap1_err, swap2_err, swap3_err, swap4_err, array_swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4
	

def teste1():
	mean_quality = 0
	mean_time = 0
	
	padrao = "/home/marilia/Dropbox/TCC - AGRUPAMENTO/Novas Instâncias - Múltiplas Habilidades/50Vertices/"
	instancia = "Categoria C/Classe8/E/IV/"
	
	S = padrao + instancia + "/S.txt"
	K = padrao + instancia + "/K.txt"
	R = padrao + instancia + "/R.txt"
	D = padrao + instancia + "/D.txt"

	data, optimal = read_data(S, K, R, D)
	
	###
	'''
	solution = {0: {0: {0.25: [], 0.5: [], 0.75: [], 1.0: [48]}, 1: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 2: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 3: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 4: {0.25: [], 0.5: [], 0.75: [], 1.0: []}}, 1: {0: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 1: {0.25: [], 0.5: [], 0.75: [2], 1.0: [8, 3, 5]}, 2: {0.25: [], 0.5: [1], 0.75: [], 1.0: []}, 3: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 4: {0.25: [], 0.5: [], 0.75: [], 1.0: []}}, 2: {0: {0.25: [0], 0.5: [], 0.75: [], 1.0: [33]}, 1: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 2: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 3: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 4: {0.25: [], 0.5: [], 0.75: [], 1.0: []}}, 3: {0: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 1: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 2: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 3: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 4: {0.25: [], 0.5: [], 0.75: [], 1.0: [42]}}, 4: {0: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 1: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 2: {0.25: [], 0.5: [], 0.75: [], 1.0: []}, 3: {0.25: [], 0.5: [], 0.75: [], 1.0: [9]}, 4: {0.25: [], 0.5: [], 0.75: [], 1.0: []}}}

	print global_efficiency(solution, S, R, 5)
	return
	'''
	###
	
	
	arq = open('/home/marilia/Dropbox/TCC - AGRUPAMENTO/Novas Instâncias - Múltiplas Habilidades/50Vertices/resultado2.txt', 'w')
	
	for i in range(0,5):
		inicio = time.time()
		
		solution, quality, swap1_err, swap2_err, swap3_err, swap4_err, swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4 = genetic_algorithm(data, optimal)
		
		
		print(solution)
		print(quality)
		return
		
		
		mean_quality += quality
		
		fim = time.time()
		
		mean_time += (fim - inicio)
	
	mean_quality /= 5
	mean_time /= 5
	
	texto = []
	texto.append('INSTANCIA = ')
	texto.append(instancia + "\n")
	texto.append('Qualidade média = ' + str(mean_quality) + "\n")
	texto.append('Tempo médio = ' + str(mean_time) + "\n")
	texto.append('\n')
	
	arq.writelines(texto)
					
	
	arq.close()


# Novas instancias, das multiplas habilidades
def multiple_skills():
	padrao = "/home/marilia/Dropbox/TCC - AGRUPAMENTO/Novas Instâncias - Múltiplas Habilidades/50Vertices/"
	categoria = ""
	classe = ""
	letra = ""
	grupo = ""
	
	S = ""
	K = ""
	R = ""
	D = ""
	
	mean_quality = 0
	mean_time = 0
	
	arq = open('/home/marilia/Dropbox/TCC - AGRUPAMENTO/Novas Instâncias - Múltiplas Habilidades/50Vertices/resultado2.txt', 'w')
	
	for c in range(0,3):
		if c == 0:
			categoria = "Categoria A"
		elif c == 1:
			categoria = "Categoria B"
		else:
			categoria = "Categoria C"
	
		for cl in range(0,3):
			if cl == 0:
				classe = "Classe2"
			elif cl == 1:
				classe = "Classe5"
			else:
				classe = "Classe8"
		
			for g in range(0,3):
				if g == 0:
					grupo = "II"
				elif g == 1:
					grupo = "III"
				else:
					grupo = "IV"
			
				for l in range(0,5):
					if l == 0:
						letra = "A"
					elif l == 1:
						letra = "B"
					elif l == 2:
						letra = "C"
					elif l == 3:
						letra = "D"
					else:
						letra = "E"
					
					mean_quality = 0
					mean_time = 0
					
					for i in range(0,5):
						inicio = time.time()
						
						S = padrao + categoria + "/" + classe + "/" + letra + "/" + grupo + "/S.txt"
						K = padrao + categoria + "/" + classe + "/" + letra + "/" + grupo + "/K.txt"
						R = padrao + categoria + "/" + classe + "/" + letra + "/" + grupo + "/R.txt"
						D = padrao + categoria + "/" + classe + "/" + letra + "/" + grupo + "/D.txt"

						data, optimal = read_data(S, K, R, D)

						solution, quality, swap1_err, swap2_err, swap3_err, swap4_err, swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4 = genetic_algorithm(data, optimal)
						
						mean_quality += quality
						
						fim = time.time()
						
						mean_time += (fim - inicio)
					
					mean_quality /= 5
					mean_time /= 5
					
					texto = []
					texto.append('INSTANCIA = ')
					texto.append(categoria + "/" + classe + "/" + letra + "/" + grupo + "\n")
					texto.append('Qualidade média = ' + str(mean_quality) + "\n")
					texto.append('Tempo médio = ' + str(mean_time) + "\n")
					texto.append('\n')
					
					arq.writelines(texto)
					
	
	arq.close()


# Instancias com uma habilidade
def one_skill():
	padrao = "/home/marilia/Dropbox/TCC - AGRUPAMENTO/instancias - CLAIO/50Vertices/"
	
	classe = ""
	letra = ""
	grupo = ""
	
	S = ""
	K = ""
	R = ""
	D = ""
	
	mean_quality = 0
	mean_time = 0
	
	arq = open('/home/marilia/Dropbox/TCC - AGRUPAMENTO/instancias - CLAIO/50Vertices/resultado-testar-swap3.txt', 'w')
	
	for cl in range(0,3):
		if cl == 0:
			classe = "Classe2"
		elif cl == 1:
			classe = "Classe5"
		else:
			classe = "Classe8"
	
		for g in range(0,3):
			if g == 0:
				grupo = "II"
			elif g == 1:
				grupo = "III"
			else:
				grupo = "IV"
		
			for l in range(0,5):
				if l == 0:
					letra = "A"
				elif l == 1:
					letra = "B"
				elif l == 2:
					letra = "C"
				elif l == 3:
					letra = "D"
				else:
					letra = "E"
				
				mean_quality = 0
				mean_time = 0
				
				print('Instancia = ' + classe + "/" + letra + "/" + grupo)
				
				for i in range(0,5):
					inicio = time.time()
					
					S = padrao + classe + "/" + letra + "/" + grupo + "/S.txt"
					K = padrao + classe + "/" + letra + "/" + grupo + "/K.txt"
					R = padrao + classe + "/" + letra + "/" + grupo + "/R.txt"
					D = padrao + classe + "/" + letra + "/" + grupo + "/D.txt"

					data, optimal = read_data(S, K, R, D)

					solution, quality, swap1_err, swap2_err, swap3_err, swap4_err, swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4 = genetic_algorithm(data, optimal)
					
					mean_quality += quality
					
					fim = time.time()
					
					mean_time += (fim - inicio)
				
				mean_quality /= 5
				mean_time /= 5
				
				texto = []
				texto.append('INSTANCIA = ')
				texto.append(classe + "/" + letra + "/" + grupo + "\n")
				texto.append('Qualidade média = ' + str(mean_quality) + "\n")
				texto.append('Tempo médio = ' + str(mean_time) + "\n")
				texto.append('\n')
				
				arq.writelines(texto)
				
				print('----------------------Instancia Finalizada----------------------------------')
	
	arq.close()

def main():
	
	#one_skill()
	#multiple_skills()
	teste1()
	
				
	#data, optimal = read_data()
	
	'''
	solution, quality, swap1_err, swap2_err, swap3_err, swap4_err, swaps_prob, qt_swap1, qt_swap2, qt_swap3, qt_swap4 = genetic_algorithm(data, optimal)
	#genetic_algorithm(data, optimal)
	
	
	if solution != None and quality != None:	
		
		print("------------------------------------------------------HORA DA VDD CARAAAAAAAAAAAI---------------------------------")
		print("QUALIDADE FINAL: " + str(quality))
		
		gap = (optimal - quality) / optimal
		print("GAP: " + str(gap))
		
		print("Probabilidade dos swaps: " + str(swaps_prob))
		print("Swap 1: Escolhido: " + str(qt_swap1) + ". Erros: " + str(swap1_err))
		print("Swap 2: Escolhido: " + str(qt_swap2) + ". Erros: " + str(swap2_err))
		print("Swap 3: Escolhido: " + str(qt_swap3) + ". Erros: " + str(swap3_err))
		print("Swap 4: Escolhido: " + str(qt_swap4) + ". Erros: " + str(swap4_err))
	
	'''
	# print("---SOLUÇÕES---")

	# solutions = initials_solutions(S, K, R, D, n, f, m, t, alfa)
	
	# swap4(solutions[0], solutions[1], m, D, n, K)
	
	'''
	for sol in solutions:
		print(sol)
		print("Qualidade: " + str(global_efficiency(sol, S, R, alfa, W)))
		
		
		new_sol = swap1(sol, R, m, f, D)
		#new_sol = swap2(sol, R, m, f, D, K)
		#new_sol = swap3(sol, alfa, R, m, f)
		
		if new_sol != None:
			print(new_sol)
			print("Qualidade: " + str(global_efficiency(new_sol, S, R, alfa, W)))
		
		
		print("-----------------------------------------------------------------------------------------------")
'''

def teste(S, R, m):
	solution = {0: {0: {0.5: [8], 1.0: [45, 49, 6, 47]}, 1: {0.5: [], 1.0: [21, 19, 25, 23]}, 2: {0.5: [], 1.0: []}, 3: {0.5: [], 1.0: []}, 4: {0.5: [], 1.0: []}}, 1: {0: {0.5: [], 1.0: [48, 20, 4, 24, 22]}, 1: {0.5: [], 1.0: []}, 2: {0.5: [], 1.0: []}, 3: {0.5: [], 1.0: []}, 4: {0.5: [], 1.0: []}}, 2: {0: {0.5: [], 1.0: [0]}, 1: {0.5: [], 1.0: [1, 3, 5, 17]}, 2: {0.5: [28], 1.0: [27]}, 3: {0.5: [], 1.0: []}, 4: {0.5: [], 1.0: []}}, 3: {0: {0.5: [], 1.0: []}, 1: {0.5: [], 1.0: []}, 2: {0.5: [], 1.0: []}, 3: {0.5: [], 1.0: []}, 4: {0.5: [], 1.0: [42, 43]}}, 4: {0: {0.5: [], 1.0: []}, 1: {0.5: [], 1.0: []}, 2: {0.5: [], 1.0: []}, 3: {0.5: [33], 1.0: [35, 34, 36, 37, 38]}, 4: {0.5: [], 1.0: []}}}
	
	a = 1 - project_efficiency(solution, 0, S, R)
	b = 1 - project_efficiency(solution, 1, S, R)
	c = 1 - project_efficiency(solution, 2, S, R)
	d = 1 - project_efficiency(solution, 3, S, R)
	e = 1 - project_efficiency(solution, 4, S, R)
	print(a + b + c + d + e)
	print(d)
	
	proj_prob = project_mutation_probability(solution, S, R, m)
	print(proj_prob)

if __name__ == "__main__":
    main()	
	
	
