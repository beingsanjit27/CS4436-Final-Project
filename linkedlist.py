# Node and Linkedlist
# ==============================
# CS4436 Final Project
# --- Sanjit Sharma --- Peter Wu --- Stefan Pisic ---
 
 
class Node:
 
    # Function to initialise the node object
    def __init__(self, dna, sequence):
        self.dna = dna  # Assign dna
        self.sequence = sequence # Assogm sequence
        self.next = None  # Initialize next as null
 
 
# Linked List class contains a Node object
class LinkedList:
 
    # Function to initialize head
    def __init__(self, head):
        self.head = head
        self.tail = head # head as tail initially
 
    # This function traverse contents of linked list
    # starting from head
    def traverseList(self):
        temp_list = []
        temp = self.head
        while (temp):
            temp_list.append(temp)
            temp = temp.next
        return temp_list


    def get_len(self):
        count = 0
        temp = self.head
        while (temp):
            count += 1
            temp = temp.next
        return count
 
 
# Code execution starts here
if __name__ == '__main__':
 
    # Start with the empty list
    llist = LinkedList()
 
    llist.head = Node(1)
    second = Node(2)
    third = Node(3)
 
    llist.head.next = second  # Link first node with second
    second.next = third  # Link second node with the third node
 
    llist.printList()