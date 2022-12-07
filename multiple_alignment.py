# Iterative multiple sequence alignment
# ==============================
# CS4436 Final Project
# --- Sanjit Sharma --- Peter Wu --- Stefan Pisic ---

from linkedlist import LinkedList, Node
from spanning_tree import Graph
import difflib

# add space to realign
def realign_sub_MPA (node_list, u_sub_MPA, v_sub_MPA, u , v, u_pairwise, v_pairwise):

    # identify which sequence is shorter
    # will add space to the shorter sequence
    # according to what "should be deleted in the longer sequence"
    u_sequence = node_list[u].sequence
    v_sequence = node_list[v].sequence

    long_sequence = u_sequence
    long_pairwise = u_pairwise
    long_sub_MPA = u_sub_MPA
    short_sequence = v_sequence
    short_pairwise = v_pairwise
    short_sub_MPA = v_sub_MPA

    if (len(u_sequence) < len(v_sequence)):
        long_sequence = v_sequence
        long_pairwise = v_pairwise
        long_sub_MPA = v_sub_MPA
        short_sequence = u_sequence
        short_pairwise = u_pairwise
        short_sub_MPA = u_sub_MPA

    # char difference between long_sequence and long_pairwise
    for i,s in enumerate(difflib.ndiff(long_sequence, long_pairwise)):
            if s[0]==' ': continue
            elif s[0]=='-': 
                
                # add the "should be deleted gap" to shorter sequence group
                for node_index in short_sub_MPA:
                    llist = list(node_list[node_index].sequence)
                    llist.insert(i, '-')
                    node_list[node_index].sequence = ''.join(llist)
                

            elif s[0]=='+':
                
                # add the "should be added gap" to longer sequence group
                for node_index in long_sub_MPA:
                    llist = list(node_list[node_index].sequence)
                    llist.insert(i, '-')
                    node_list[node_index].sequence = ''.join(llist)


    # char difference between short_sequence and short_pairwise
    for i,s in enumerate(difflib.ndiff(short_sequence, short_pairwise)):
            if s[0]==' ': continue
            elif s[0]=='-':    # add the "should be deleted gap" to longer sequence group
                for node_index in long_sub_MPA:
                    llist = list(node_list[node_index].sequence)
                    llist.insert(i, '-')
                    node_list[node_index].sequence = ''.join(llist)

            elif s[0]=='+':
                
                # add the "should be added gap" to shorter sequence group
                for node_index in short_sub_MPA:
                    llist = list(node_list[node_index].sequence)
                    llist.insert(i, '-')
                    node_list[node_index].sequence = ''.join(llist)

    # for node in node_list:
    #     print(node.dna)
    #     print(node.sequence)

    return node_list

def iterative_alignment(tree_order, pairwise_alignments, dna_dict_num):

    # create node for each dna
    # list to store nodes
    node_list = []
    MPA = [] # Multiple Process Alignment
    for dna in dna_dict_num.keys():
        node = Node(dna, '')
        node_list.append(node)

    # when adding new node to list, add from tail, and make new node as tail
    # also remove the old linkedlist
    for u, v in tree_order:

        u_dna = [i for i in dna_dict_num if dna_dict_num[i]==u][0]
        v_dna = [i for i in dna_dict_num if dna_dict_num[i]==v][0]

        # check if either node is not in MPA 
        u_in = False
        v_in = False
        u_sub_MPA = []
        v_sub_MPA = []
        if len(MPA) > 0:
            for i in MPA:
                if u in i:
                    u_in = True
                    u_sub_MPA = i
                if v in i:
                    v_in = True
                    v_sub_MPA = i

        # both not in MPA
        if ( not u_in and not v_in):

            for alignment in pairwise_alignments:

                if (alignment[0] == u_dna and alignment[1] == v_dna):
                    node_list[u].sequence = alignment[3]
                    node_list[v].sequence = alignment[4]
                    
                elif (alignment[1] == u_dna and alignment[0] == v_dna):
                    node_list[u].sequence = alignment[4]
                    node_list[v].sequence = alignment[3]

             # add pairwise as a new multiple alignment group
            MPA.append([u,v]) 

        # u in MPA but v not in
        # add v from pairwise
        elif (u_in and not v_in):
            for alignment in pairwise_alignments:
                if (alignment[0] == u_dna and alignment[1] == v_dna):

                    node_list[v].sequence = alignment[4] # get v's sequecne
                    # realign v and u's group
                    node_list = realign_sub_MPA(node_list, u_sub_MPA, [v], u, v, alignment[3], alignment[4])
                    
                elif (alignment[1] == u_dna and alignment[0] == v_dna):

                    node_list[v].sequence = alignment[3] # get v's sequecne
                    # realign v and u's group
                    node_list = realign_sub_MPA(node_list, u_sub_MPA, [v], u, v, alignment[4], alignment[3])
            

            # remove old sub_MPA
            MPA.remove(u_sub_MPA)
            # add new sub_MPA
            u_sub_MPA.append(v)
            MPA.append(u_sub_MPA)

        # u not in MPA but v in
        elif (not u_in and v_in):
            for alignment in pairwise_alignments:
                if (alignment[0] == u_dna and alignment[1] == v_dna):

                    node_list[u].sequence = alignment[3] # get u's sequecne
                    # realign u and v's group
                    node_list = realign_sub_MPA(node_list, v_sub_MPA, [u], v, u, alignment[4], alignment[3])
                    
                elif (alignment[1] == u_dna and alignment[0] == v_dna):

                    node_list[u].sequence = alignment[4] # get u's sequecne
                    # realign u and v's group
                    node_list = realign_sub_MPA(node_list, v_sub_MPA, [u], v, u, alignment[3], alignment[4])
            

            # remove old sub_MPA
            MPA.remove(v_sub_MPA)
            # add new sub_MPA
            v_sub_MPA.append(u)
            MPA.append(v_sub_MPA)

        # both already in MPA
        else:
            for alignment in pairwise_alignments:
                if (alignment[0] == u_dna and alignment[1] == v_dna):
                     # realign u's and v's group
                    node_list = realign_sub_MPA(node_list, u_sub_MPA, v_sub_MPA, u, v, alignment[3], alignment[4])
                elif (alignment[1] == u_dna and alignment[0] == v_dna):
                    # realign u's and v's group
                    node_list = realign_sub_MPA(node_list, u_sub_MPA, v_sub_MPA, u, v, alignment[4], alignment[3])

            # adjust MPA
            MPA.remove(u_sub_MPA)
            MPA.remove(v_sub_MPA)
            MPA.append( v_sub_MPA + u_sub_MPA)
    
    for node in node_list:
        print(node.dna)

    for node in node_list:
        #print(node.dna)
        print(node.sequence.replace("\n", ""))

    # pairwise sequence alignment
    # store alignment results
    with open('multiple_results.txt', 'w') as f:
        for node in node_list:

            f.write(node.dna + '\n')
            f.write(node.sequence)

    return