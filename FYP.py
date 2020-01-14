import tkinter as tk
from tkinter import *
from tkinter import Tk
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
# import matplotlib.pyplot as plt


mw = tk.Tk()  # type:Tk
var = IntVar()
def openfile():
    global my_seq
    global antigenicity
    global m, a, c, b
    from tkinter import filedialog
    root = Tk()
    root.filename = filedialog.askopenfilename(initialdir="/", title="Select file",
                                               filetypes=(("pdb files", "*.pdb"), ("pdb files", "*.pdb")))
    print(root.filename)
    structure_id = "1e6j"
    structure = PDBParser().get_structure(structure_id, root.filename)
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure):
        my_seq = pp.get_sequence()  # type: Seq
        print(my_seq)
    for model in structure:
        for chain in model:
            print(chain)
    sequence = list(my_seq)
    m = ''.join(sequence)  # type: str
    print(m)
    length = len(m)  # type: int
    print(length)
    print("Sequence consist of", len(m), "Amino Acids")
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    analysed_seq = ProteinAnalysis(m)
    print("Molecular weight = ", analysed_seq.molecular_weight())
    print("Amino Acid Count = ", analysed_seq.count_amino_acids())
    print("Secondary structure fraction =", analysed_seq.secondary_structure_fraction())
    kd = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
          'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
          'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
          'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    c = list(analysed_seq.flexibility())
    b = list(analysed_seq.protein_scale(kd, 10, 1.0))
    i = 1
    j = -1  # type: int
    k = 9
    while i <= (length - 10):
        print("Sequence is = ", m[j + 1:k + 1])
        print("Flexibility value = ", c[j + 1])
        print("Hydrophilicity value = ", b[j + 1])
        ana_seq = ''.join(m[j + 1:k + 1])
        analyze_seq = ProteinAnalysis(ana_seq)
        # For Secondary structure Analysis
        print("Secondary structure fraction =", analyze_seq.secondary_structure_fraction())
        a = list(analyze_seq.secondary_structure_fraction())
        global tupleall
        tupleall = (m[j + 1:k + 1], c[j + 1], b[j + 1], a)
        print(tupleall[0], tupleall[2], tupleall[1], tupleall[3])
        i = i + 1
        if a[0] >= a[1]:
            a[0] = 1
        else:
            a[0] = a[1]
        # For Hydrophilicity
        if b[j + 1] > 0.5:
            b[j + 1] = 2
        elif b[j + 1] < 0.5 or b[j + 1] > 0:
            b[j + 1] = 1
        elif b[j + 1] > 0 or b[j + 1] > -0.4:
            b[j + 1] = -1
        elif b[j + 1] < -0.4:
            b[j + 1] = -2
        else:
            b[j + 1] = 0
        # For Flexibility
        if c[j + 1] > 1.0:
            c[j + 1] = 1
        else:
            c[j + 1] = 0
        # For antigenicity Index
        antigenicity = 0.3 * b[j + 1] + 0.15 * 1 + 0.15 * c[j + 1] + 0.2 * a[0]
        print("antigenicity", antigenicity)
        j += 1
        k += 1


def result():
    global item
    item = str(tupleall[0]) + str(tupleall[2]) + str(tupleall[1]) + str(tupleall[3])
    # if you want the button to disappear:
    # button.destroy() or button.pack_forget()
    label = Label(mw, text=item, bg="black", fg="white")
    # this creates a new label to the GUI
    label.pack()


mw.option_add("*Button.Background", "black")
mw.option_add("*Button.Foreground", "white")
mw.title('Protein Antigenicity Finder')
mw.geometry("1300x600")  # You want the size of the app to be 1300x600
mw.resizable(0, 0)  # Don't allow resizing in the x or y direction
back = tk.Frame(master=mw, bg='black')
back.pack_propagate(0)  # Don't allow the widgets inside to determine the frame's width / height
back.pack(fill=tk.BOTH, expand=1)  # Expand the frame to fill the root window
e = tk.Label(master=back, text='Upload PDB File', bg="black", fg="white",
             font=('Times New Roman', 18, 'bold')).grid(row=1, column=1)
n = tk.Button(master=back, text='Upload File', bg="grey", fg="black", command=openfile).grid(row=3, column=1)
info0 = tk.Label(master=back, text='Functionalities', bg="black", fg="green",
                 font=('Times New Roman', 14, 'bold')).grid(row=6, column=1)
info0 = tk.Label(master=back, text='Hydrophilicity', bg="black", fg="pink",
                 font=('Times New Roman', 14, 'bold')).grid(row=7, column=1)
info0 = tk.Label(master=back, text='Flexibility', bg="black", fg="pink",
                 font=('Times New Roman', 14, 'bold')).grid(row=8, column=1)
info0 = tk.Label(master=back, text='Accessibility', bg="black", fg="pink",
                 font=('Times New Roman', 14, 'bold')).grid(row=8, column=2)
info0 = tk.Label(master=back, text='Secondary Structure Value', bg="black", fg="pink",
                 font=('Times New Roman', 14, 'bold')).grid(row=7, column=2)
info9 = tk.Button(master=back, text='Download Result', bg="green", fg="white", command=result).grid(row=12, column=1)
info10 = tk.Label(master=back, text='', bg="black", fg="white").grid(row=13, column=1)
info01 = tk.Label(master=back, text="Sequence   Hydrophlicity    Flexibility     Secondary structure ", bg="black",
                  fg="white",
                  font=('Times New Roman', 15, 'bold')).grid(row=15, column=2)
info = tk.Label(master=back, bg="black")
close = tk.Button(master=back, text='Quit', bg="red", fg="white", command=mw.destroy).grid(row=12, column=2)
info.pack()
mw.mainloop()
