# program main
# ==============================
# CS4436 Final Project
# --- Sanjit Sharma --- Peter Wu --- Stefan Pisic ---

from itertools import combinations

from multiple_alignment import iterative_alignment
from sequence_alignment import nw, score_matrix
from spanning_tree import Graph


# read alignments
# alignments[x][0]: dna1
# alignments[x][1]: dna2
# alignments[x][2]: alignment score
# alignments[x][3]: alignment for sequence 1 
# alignments[x][4]: alignment for sequence 2
def read_alignments(file):
    alignments = []
    with open('alignment_results.txt', 'r') as f:
        dnas = f.readline()
        while(True):
            dnas = dnas.split()
            score = f.readline().split()
            sequence_1 = f.readline()
            sequence_2 = f.readline()
            alignments.append([dnas[0], dnas[1], float(score[0]), sequence_1, sequence_2])

            dnas = f.readline()
            if dnas == '':
                break
            
    return alignments

# main
if __name__ == '__main__':

    human = 'VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'
    olive_baboon = 'GLSDGEWQLVLNVWGKVEADIPSHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLELISESIIQVLQSKHPGDFGADAQGAMNKALELFRNDMAAKYKELGFQG'
    garden_pea = 'GFTDKQEALVNSSSEFKQNLPGYSILFYTIVLEKAPAAKGLFSFLKDTAGVEDSPKLQAHAEQVFGLVRDSAAQLRTKGEVVLGNATLGAIHVQKGVTNPHFVVVKEALLQTIKKASGNNWSEELNTAWEVAYDGLATAIKKAMKTA'
    gayal = 'VLSAADKGNVKAAWGKVGDHAAEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGAKVAAALTKAVGHLDDLPGALSELSDLHAHKLRVDPVNFKLLSHSLLVTLASHLPNDFTPAVHASLDKFLANVSTVLTSKYR'
    goat  = 'VLSAADKSNVKAAWGKVGGNAGAYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGEKVAAALTKAVGHLDDLPGTLSDLSDLHAHKLRVDPVNFKLLSHSLLVTLACHLPNDFTPAVHASLDKFLANVSTVLTSKYR'
    mouse = 'VHLTDAEKSAVSCLWAKVNPDEVGGEALGRLLVVYPWTQRYFDSFGDLSSASAIMGNPKVKAHGKKVITAFNEGLKNLDNLKGTFASLSELHCDKLHVDPENFRLLGNAIVIVLGHHLGKDFTPAAQAAFQKVVAGVATALAHKYH'
    sheep = 'MLTAEEKASVISLFAKVNVEEVGGEALGRLLVVYPWTQRFFEHFGDLSSADAILGNPKVKGHGKKVLNSFSEGLKQLDDLKGAFASLSELHCDKLHVDPENFRLLGNVLVVVLARRFGGEFTPELQANFQKVVTGVANALAHRYH'
    white_rhinoceres = 'VELTAEEKAAVLALWDKVKEDEVGGEALGRLLVVYPWTQRFFDSFGDLSTPAAVMGNAKVKAHGKKVLHSFGDGVHHLDNLKGTFAALSELHCDKLHVDPENFRLLGNVLVVVLAKHFGKQFTPELQAAYQKVVAGVANALAHKYH'
    european_moose = 'VLSATDKSNVKAAWGKVGGNAPAYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKAHGEKVANALTKAVGHLDDLPGTLSDLSDLHAHKLRVDPVNFKLLSHTLLVTLAAHLPSDFTPAVHASLDKFLANVSTVLTSKYR'
    sesbania_rostrata = 'MGFTEKQEALVNASYEAFKQNLPGNSVLFYSFILEKAPAAKGMFSFLKDSDGVPQNNPSLQAHAEKVFGLVRDSAAQLRATGVVVLADASLGSVHVQKGVLDPHFVVVKEALLKTLKEAAGATWSDEVSNAWEVAYDGLSAAIKKAMS'

    dna_dict= {
        'human': human,
        'olive_baboon': olive_baboon,
        'garden_pea': garden_pea,
        'gayal': gayal,
        'goat': goat,
        'mouse': mouse,
        'sheep': sheep,
        'white_rhinoceres': white_rhinoceres,
        'european_moose': european_moose,
        'sesbania_rostrata': sesbania_rostrata
    }
    dna_pairs =  combinations(dna_dict.keys(), 2) # non-orderd pairs of DNAs

    blosum62_dict = score_matrix('CSTPAGNDEQHRKMILVFYW')    

    # pairwise sequence alignment
    # store alignment results
    with open('alignment_results.txt', 'w') as f:
        for dna1, dna2 in dna_pairs:

            # Needleman-Wunsch Algorithm/ Sequence Alignment
            alignments = nw(dna_dict[dna1], dna_dict[dna2], blosum62_dict)

            score  = alignments[0]
            alignment = alignments[1]

            f.write(dna1 + ' ' + dna2+'\n')
            f.write(str(score) + '\n')
            f.write(alignment + '\n')

    print('####pairwise sequence alignment done')

    pairwise_alignments = read_alignments('alignment_results.txt')

    dnas = ['human', 'olive_baboon', 'garden_pea', 'gayal', 'goat', 'mouse', 'sheep',
    'white_rhinoceres', 'european_moose', 'sesbania_rostrata']

    # encode dnas into graph node
    g = Graph(len(dnas)) 
    nodes = range(0, len(dnas))
    dna_dict_num = dict()
    for i in nodes:
        dna_dict_num[dnas[i]] = i

    for i in pairwise_alignments:

        u = dna_dict_num.get(i[0])
        v = dna_dict_num.get(i[1])
        g.addEdge(u, v, -i[2]) # negative score for max spanning tree
 
    # Kruska'ss min spanning tree
    result = g.KruskalMST()
    tree_order = []
    for u,v, weight in result:
        tree_order.append([u, v])

    print('####spanning tree to guide the order of multiple alignment done')

    # Iterative multiple alignment
    iterative_alignment(tree_order, pairwise_alignments, dna_dict_num)

    print('####iterative multiple alignment done')

    import difflib