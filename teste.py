import os
from collections import defaultdict
import matplotlib.pyplot as plt
import random

# Função para carregar sequência de um arquivo FASTA
def load_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        sequence = ''.join(lines[1:]).replace('\n', '')
    return sequence

# Função para contar os aminoácidos e suas posições
def analise_sequencia(sequencia):
    aminoacidos = defaultdict(lambda: {'count': 0, 'positions': []})
    for i, aa in enumerate(sequencia):
        aminoacidos[aa]['count'] += 1
        aminoacidos[aa]['positions'].append(i+1)
    return aminoacidos

# Função para escrever a análise em um arquivo
def write_analysis(filename, aminoacidos):
    with open(filename, 'w') as file:
        file.write("Código do aminoácido\tNome completo\tNúmero de ocorrências\tPosições\n")
        for aa, data in aminoacidos.items():
            file.write(f"{aa}\t{aminoacido_nomes[aa]}\t{data['count']}\t{','.join(map(str, data['positions']))}\n")

# Dicionário com os nomes completos dos aminoácidos
aminoacido_nomes = {
    'A': 'Alanina', 'R': 'Arginina', 'N': 'Asparagina', 'D': 'Ácido aspártico', 'C': 'Cisteína',
    'E': 'Ácido glutâmico', 'Q': 'Glutamina', 'G': 'Glicina', 'H': 'Histidina', 'I': 'Isoleucina',
    'L': 'Leucina', 'K': 'Lisina', 'M': 'Metionina', 'F': 'Fenilalanina', 'P': 'Prolina',
    'S': 'Serina', 'T': 'Treonina', 'W': 'Triptofano', 'Y': 'Tirosina', 'V': 'Valina'
}

# Carregar e analisar sequências de diferentes organismos
organismos = ['organism1.fasta', 'organism2.fasta', 'organism3.fasta', 'organism4.fasta']
for organismo in organismos:
    sequencia = load_fasta(organismo)
    aminoacidos = analise_sequencia(sequencia)
    write_analysis(f'{os.path.splitext(organismo)[0]}-aminoacidos.txt', aminoacidos)

# Criar histogramas de ocorrências de aminoácidos
for organismo in organismos:
    sequencia = load_fasta(organismo)
    aminoacidos = analise_sequencia(sequencia)
    counts = [data['count'] for data in aminoacidos.values()]
    plt.figure(figsize=(10, 6))
    plt.bar(aminoacidos.keys(), counts)
    plt.xlabel('Aminoácidos')
    plt.ylabel('Número de Ocorrências')
    plt.title(f'Número de Ocorrências de Aminoácidos em {os.path.splitext(organismo)[0]}')
    plt.savefig(f'{os.path.splitext(organismo)[0]}-histograma.png')
    plt.close()

def compare_sequences(seq1, seq2):
    differences = []
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            differences.append((i+1, seq1[i], seq2[i]))
    return differences

# Carregar e comparar sequências
sequences = [load_fasta(f'organism{i+1}.fasta') for i in range(4)]
comparisons = [(1, 2), (1, 3), (1, 4)]  # Comparações entre pares de organismos

for org1, org2 in comparisons:
    differences = compare_sequences(sequences[org1-1], sequences[org2-1])
    with open(f'comparacao_organism{org1}_vs_organism{org2}.txt', 'w') as file:
        file.write("Posição\tOrganismo 1\tOrganismo 2\n")
        for diff in differences:
            file.write(f"{diff[0]}\t{diff[1]}\t{diff[2]}\n")

# Função para gerar uma mutação pontual na sequência
def mutate_sequence(sequence):
    position = random.randint(0, len(sequence) - 1)
    new_aa = random.choice(list(aminoacido_nomes.keys()))
    mutated_sequence = sequence[:position] + new_aa + sequence[position + 1:]
    return mutated_sequence, position + 1, new_aa

# Carregar a sequência original e gerar uma mutação
for i in range(4):
    sequence = load_fasta(f'organism{i+1}.fasta')
    mutated_sequence, position, new_aa = mutate_sequence(sequence)
    with open(f'organism{i+1}-mutacao.fasta', 'w') as file:
        file.write(f">organism{i+1}-mutacao\n")
        file.write('\n'.join([mutated_sequence[j:j+60] for j in range(0, len(mutated_sequence), 60)]))

    differences = compare_sequences(sequence, mutated_sequence)
    with open(f'comparacao_organism{i+1}_vs_mutacao.txt', 'w') as file:
        file.write("Posição\tOrganismo 1\tMutacao\n")
        for diff in differences:
            file.write(f"{diff[0]}\t{diff[1]}\t{diff[2]}\n")
