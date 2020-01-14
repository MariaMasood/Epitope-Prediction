import tkinter as tk
from tkinter import *
from tkinter import Tk
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Seq import Seq
import freesasa
import math
mw = tk.Tk()  # type:Tk
var = IntVar()


def openfile():
    global prob, probab, te
    global my_seq
    global anti
    global structure, structure_id, filename
    global antigenicity, hydro, flex, sec
    global m, a, c, b, length, j, k
    global hydroph, flexi, access
    anti = []
    sec = []
    probab = []
    from tkinter import filedialog
    root = Tk()
    root.filename = filedialog.askopenfilename(initialdir="/", title="Select file",
                                               filetypes=(("pdb files", "*.pdb"), ("pdb files", "*.pdb")))
    filename = root.filename
    print(filename)
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
    m = ''.join(sequence)
    print(m)
    length = len(m)  # type: int
    print("Sequence consist of", length, "Amino Acids")
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
    hydro = list(analysed_seq.protein_scale(kd, 10, 1.0))
    flex = list(analysed_seq.flexibility())
    hydroph = list(analysed_seq.protein_scale(kd, 10, 1.0))
    flexi = list(analysed_seq.flexibility())

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
        a = a[0]
        sec.append(a)
        i += 1
        j += 1
        k += 1
    f = length
    r = 1
    y = 10
    global acc, logacc
    acc = []
    for i in range(0, f):
        str1 = "accessibility, resi "
        str2 = str(r) + "-" + str(y)
        saving = str1 + str2
        print(saving)
        r = r + 1
        y = y + 1
        structure = freesasa.Structure("1e6j.pdb")
        resulta = freesasa.calc(structure)
        area_classes = freesasa.classifyResults(resulta, structure)
        print("Total : %.2f A2" % resulta.totalArea())
        for key in area_classes:
            print(key, ": %.2f A2" % area_classes[key])
        resulta = freesasa.calc(structure, freesasa.Parameters({'algorithm': freesasa.LeeRichards, 'n-slices': 10}))
        selections = freesasa.selectArea(('alanine, resn ala', saving), structure, resulta)
        for key in selections:
            print(key, ": %.2f A2" % selections[key])
            a = selections[key]
            acc.append(a)

    l = acc[0::2]
    access = l
    print(acc)
    print(l)
    logacc = [math.log(y, 10) for y in l]

    print(logacc)


def result():
    global item
    item = ["Your result downloaded in Result folder"]
    download_dir = "Result\Antigenicity.csv"  # where you want the file to be downloaded to
    csv = open(download_dir, "w")
    # "w" indicates that you're writing strings to the file
    columnTitleRow = 'Sequence, Hydrophilicity, Flexibility, Log Accessibility, Antigenicity,Probable Antigen \n '
    csv.write(columnTitleRow)
    # for key in dic.keys():
    sequence = list(my_seq)
    m = ''.join(sequence)  # type: str
    i = 1
    while i <= (130):
        i += 1
       # if a[0] >= a[1]:
          #  a[0] = 1
       # else:
           # a[0] = a[1]
        # For Hydrophilicity
        if hydro[i] > 0.5:
            hydro[i] = 2
        elif hydro[i] < 0.5 or hydro[i] > 0:
            hydro[i] = 1
        elif hydro[i] > 0 or hydro[i] > -0.4:
            hydro[i] = -1
        elif hydro[i] < -0.4:
            hydro[i] = -2
        else:
            hydro[i] = 0
        # For Flexibility
        if flex[i] > 1.0:
            flex[i] = 1
        else:
            flex[i] = 0
        # For accessibility
        if logacc[i] > 1.0:
            logacc[i] = 1
        else:
            logacc[i] = 0
            # For plrobable antigenecity

        # For antigenicity Index
        antigenicity = 0.3 * hydro[i] + 0.15 * logacc[i] + 0.15 * flex[i] + 0.2 * 1
        print("antigenicity", antigenicity)
        temp = antigenicity
        anti.append(temp)
        if temp > 1.0:
            prob = "Yes"
            te = prob
            probab.append(te)
        else:
            prob = "No"
            te = prob
            probab.append(te)





    x = map(lambda hydroph: str(hydroph), hydroph)
    x = list(x)
    print("hydro",x)
    y = map(lambda flexi: str(flexi), flexi)
    y = list(y)
    print("flex",y)
    w = map(lambda access: str(access), access)
    w = list(w)
    j = map(lambda anti: str(anti), anti)
    j = list(j)
    print(anti)
    g = -1  # type: int
    h = 9
    for i in range(130):
        row = m[g + 1:h + 1] + "," + x[i] + "," + y[i] + "," + w[i] + "," + j[i] + "," + probab[i] +"," + "\n"
        g += 1
        h += 1
        csv.write(row)
    # if you want the button to disappear:
    # button.destroy() or button.pack_forget()
    label = Label(mw, text=item, bg="black", fg="white")
    # this creates a new label to the GUI
    label.pack()


mw.option_add("*Button.Background", "grey")
mw.option_add("*Button.Foreground", "white")
mw.title('Protein Antigenicity Finder')
mw.geometry("900x600")  # You want the size of the app to be 1300x600
mw.resizable(0, 0)  # Don't allow resizing in the x or y direction
back = tk.Frame(master=mw, bg='grey')
back.pack_propagate(0)  # Don't allow the widgets inside to determine the frame's width / height
back.pack(fill=tk.BOTH, expand=1)  # Expand the frame to fill the root window
e = tk.Label(master=back, text='Baqai Medical University', bg="grey", fg="white",
             font=('Times New Roman', 18, 'bold')).grid(row=1, column=2)
e01 = tk.Label(master=back, text='Baqai Institute of Information Technology', bg="grey", fg="white",
               font=('Times New Roman', 18, 'bold')).grid(row=2, column=2)
e02 = tk.Label(master=back, text='Upload PDB File', bg="grey", fg="black",
               font=('Times New Roman', 18, 'bold')).grid(row=5, column=1)
e03 = tk.Label(master=back, text='', bg="grey", fg="white",
               font=('Times New Roman', 18, 'bold')).grid(row=4, column=3)
e04 = tk.Label(master=back, text='', bg="grey", fg="white",
               font=('Times New Roman', 18, 'bold')).grid(row=8, column=3)
n = tk.Button(master=back, text='Upload File', bg="grey", fg="black", command=openfile).grid(row=7, column=1)
info9 = tk.Button(master=back, text='Download Result', bg="green", fg="white", command=result).grid(row=12, column=1)
info10 = tk.Label(master=back, text='', bg="grey", fg="white").grid(row=13, column=1)
info = tk.Label(master=back, bg="grey")
close = tk.Button(master=back, text='Quit', bg="red", fg="white", command=mw.destroy).grid(row=12, column=2)
info.pack()
mw.mainloop()
