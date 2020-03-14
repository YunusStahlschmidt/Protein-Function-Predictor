from Tkinter import *
import tkFileDialog
from recommendations import *


class ProteinFunctionTool(Frame):
    def __init__(self, parent):  # Initializing the UI
        Frame.__init__(self)
        self.initUI()
        self.data_center = DataCenter()  # Creating the data center where all data is kept

    def initUI(self):
        # Title Label for header
        self.label_title = Label(text="Protein Function Prediction Tool", font=("", "20", "bold"), fg="white", bg="mediumspringgreen", height=2)
        self.label_title.pack(fill=X)

        # Frame that holds the buttons for
        self.frame_upload_buttons = Frame()  # Frame for the upload buttons
        self.frame_upload_buttons.pack(pady=30)
        self.button_upload_annotations = Button(self.frame_upload_buttons, text="Upload"+'\n'+"Annotations", command=self.OpenAnnoFile)
        self.button_upload_annotations.grid(row=0, column=0, padx=30)
        self.button_upload_evidence_code = Button(self.frame_upload_buttons, text="Upload Evidence"+'\n'+"Code Values", command=self.OpenECVFile)
        self.button_upload_evidence_code.grid(row=0, column=1, padx=30)
        self.button_upload_go_file = Button(self.frame_upload_buttons, text="Upload"+'\n'+"GO File", command=self.OpenGOFile)
        self.button_upload_go_file.grid(row=0, column=2, padx=30)
        # Frame for the list boxes and the checkbuttons
        self.frame_list_boxes = Frame()
        self.frame_list_boxes.pack()
        # Frame for the first listbox and scrollbar
        self.label_proteins = Label(self.frame_list_boxes, text="Proteins")
        self.label_proteins.grid(row=0, column=0)
        self.frame_proteins = Frame(self.frame_list_boxes)
        self.frame_proteins.grid(row=1, column=0)
        self.scrollbar_proteins = Scrollbar(self.frame_proteins)
        self.scrollbar_proteins.pack(side=RIGHT, fill=Y)
        self.listbox_proteins = Listbox(self.frame_proteins, yscrollcommand=self.scrollbar_proteins.set, width=28, height=15)
        self.listbox_proteins.bind("<<ListboxSelect>>", self.on_select)
        self.listbox_proteins.pack(side=LEFT)
        self.scrollbar_proteins.configure(command=self.listbox_proteins.yview)
        # Frame for the checkbuttons for setting the similarity type
        self.label_similarity_metric = Label(self.frame_list_boxes, text="Similarity Metric", padx=50)
        self.label_similarity_metric.grid(row=0, column=1)
        self.frame_similarity_metric = Frame(self.frame_list_boxes, borderwidth=5, highlightbackground="black", highlightthickness=3)
        self.frame_similarity_metric.grid(row=1, column=1)
        self.sim_var = StringVar()
        self.sim_var.set(" ")
        self.checkbutton_pearson = Checkbutton(self.frame_similarity_metric, text="Pearson", variable=self.sim_var, onvalue="Pearson", command=self.calculations)
        self.checkbutton_pearson.pack(pady=(30, 0))
        self.checkbutton_euclidean = Checkbutton(self.frame_similarity_metric, text="Euclidean", variable=self.sim_var, onvalue="Euclidean", command=self.calculations)
        self.checkbutton_euclidean.select()
        self.checkbutton_euclidean.pack(pady=(0, 210), padx=10)
        # Frame for the second listbox and scrollbar
        self.label_similar_protein = Label(self.frame_list_boxes, text="Similar Protein")
        self.label_similar_protein.grid(row=0, column=2)
        self.frame_similar_protein = Frame(self.frame_list_boxes)
        self.frame_similar_protein.grid(row=1, column=2)
        self.scrollbar_similar_protein = Scrollbar(self.frame_similar_protein)
        self.scrollbar_similar_protein.pack(side=RIGHT, fill=Y)
        self.listbox_similar_protein = Listbox(self.frame_similar_protein, yscrollcommand=self.scrollbar_similar_protein.set, width=40, height=15)
        self.listbox_similar_protein.pack(side=LEFT)
        self.scrollbar_similar_protein.configure(command=self.listbox_similar_protein.yview)
        # Frame for the third listbox and scrollbar
        self.label_predicted_function = Label(self.frame_list_boxes, text="Predicted Function", padx=40)
        self.label_predicted_function.grid(row=0, column=3)
        self.frame_predicted_function = Frame(self.frame_list_boxes, padx=40)
        self.frame_predicted_function.grid(row=1, column=3)
        self.scrollbar_predicted_function = Scrollbar(self.frame_predicted_function)
        self.scrollbar_predicted_function.pack(side=RIGHT, fill=Y)
        self.listbox_predicted_function = Listbox(self.frame_predicted_function, yscrollcommand=self.scrollbar_predicted_function.set, width=80, height=15)
        self.listbox_predicted_function.pack(side=LEFT)
        self.scrollbar_predicted_function.configure(command=self.listbox_predicted_function.yview)
        self.pack(fill=BOTH)

    def on_select(self, event):  # Calls a function to calculate the similar proteins and the predicted functions for the selected protein
        widget = self.listbox_proteins  # Gets the selected protein name from the 1. listbox
        selection = widget.curselection()
        self.prot_lb_selection = widget.get(selection)
        self.calculations()

    def calculations(self):  # Calculates the similar proteins and predicted functions and displays them if item is selected on the listbox or the similarity metric is changeds
        self.listbox_similar_protein.delete(0, END)  # Deleting the stuff from the listboxes if there is something
        self.listbox_predicted_function.delete(0, END)
        sim = self.sim_var.get()  # Gets the similarity that is set by the checkbuttons
        if sim == "Pearson":  # Setting the similarity for the functions of recommendations.py
            sim = sim_pearson
        elif sim == "Euclidean":
            sim = sim_distance
        self.similar_proteins = topMatches(self.data_center.recommendations_dict, self.prot_lb_selection, similarity=sim, n=len(self.data_center.proteins_dict))  # Calculating the similar proteins
        self.similar_proteins_rec_dict = {}  # creating a dictionary for the calculation of predicted functions                # for the selected protein by using topMatches from recommendations.py
        for item in self.similar_proteins:
            if item[0] > 0:  # Do if the similarity is bigger than 0
                insert_list = []  # List for keeping the elements of 1 row that are going to be inserted into the 2. listbox
                insert_list.append(str(item[0]))
                keys_prot_dict = self.data_center.proteins_dict.keys()
                for key in keys_prot_dict:
                    if self.data_center.proteins_dict[key].name == item[1]:  # looking up which protein is beeing handled
                        my_protein = self.data_center.proteins_dict[key]
                        insert_list.append(my_protein.id)
                        protein = self.similar_proteins_rec_dict[item[1]] = {}
                        keys_annot_dict = self.data_center.proteins_dict[key].annotations_dict.keys()
                        for anno_key in keys_annot_dict:  # This part looks for the functionality of the protein that is selected
                            self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.term_name = self.data_center.go_dict[anno_key]
                            funct = self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.term_name
                            protein[funct] = float(self.data_center.proteins_dict[key].annotations_dict[anno_key].evidence_code.numeric_value)
                        break
                insert_list.append(item[1])
                self.listbox_similar_protein.insert(END, " - ".join(insert_list))  # Inserting the similarity score, protein id and protein name to 2. listbox
            elif item[0] <= 0:  # If similerity is <= 0 breaks for loop for shorter runtime and better memory usage
                break
        self.similar_proteins_rec_dict[self.prot_lb_selection] = protein = {}  # Adding the selected protein to the dictionary to use in getRecommendations
        keys_prot_dict = self.data_center.proteins_dict.keys()
        for key in keys_prot_dict:
            if self.data_center.proteins_dict[key].name == self.prot_lb_selection:  # Finding which protein we are dealing with
                my_protein = self.data_center.proteins_dict[key]
                keys_annot_dict = self.data_center.proteins_dict[key].annotations_dict.keys()
                for anno_key in keys_annot_dict:  # Setting the proteins properties
                    self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.term_name = self.data_center.go_dict[anno_key]
                    funct = self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.term_name
                    protein[funct] = float(self.data_center.proteins_dict[key].annotations_dict[anno_key].evidence_code.numeric_value)
                break

        self.predicted_functions = getRecommendations(self.similar_proteins_rec_dict, self.prot_lb_selection, similarity=sim)  # Using getRecommendations from recommendations.py to calculate
        for item in self.predicted_functions:                                                                                  # the predicted functions of the similar proteins
            insert_list = []  # Creating a list to keep the stuff thats going to be inserted
            insert_list.append(str(item[0]))
            keys_prot_dict = self.data_center.proteins_dict.keys()
            for key in keys_prot_dict:
                keys_annot_dict = self.data_center.proteins_dict[key].annotations_dict.keys()
                for anno_key in keys_annot_dict:
                    if self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.term_name == item[1]:
                        if self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.id not in insert_list:  # Since some GO id are shared by multiple proteins we check if the
                            insert_list.append(self.data_center.proteins_dict[key].annotations_dict[anno_key].functionality.id)  # GO id already exists in the insert list
                        break
            insert_list.append(item[1])
            self.listbox_predicted_function.insert(END, " - ".join(insert_list))  # inserting the data to the 3. listbox

    def OpenAnnoFile(self):  # Calls a method from DataCenter class to open file browser to select GEO_human.gaf file
        self.data_center.RnP_Annotation()
        keys = self.data_center.proteins_dict.keys()
        for key in keys:  # Inserting all the protein names to the 1. listbox
            self.listbox_proteins.insert(END, self.data_center.proteins_dict[key].name)
        self.update_idletasks()

    def OpenECVFile(self):  # Calls a method from DataCenter class to open file browser to select ecv.txt file
        self.data_center.RnP_ECV()
        keys_prot_dict = self.data_center.proteins_dict.keys()
        for prot_key in keys_prot_dict:  # Updating the already existing protein objects ecv to its actual value
            keys_annot_dict = self.data_center.proteins_dict[prot_key].annotations_dict.keys()
            for anno_key in keys_annot_dict:
                ec_acronym = self.data_center.proteins_dict[prot_key].annotations_dict[anno_key].evidence_code.acronym
                ec_num_value = self.data_center.evidence_codes_dict[ec_acronym].numeric_value
                self.data_center.proteins_dict[prot_key].annotations_dict[anno_key].evidence_code.numeric_value = ec_num_value

    def OpenGOFile(self):  # Calls a method from DataCenter class to open file browser to select go.obo file
        self.data_center.RnP_GO()
        keys_prot_dict = self.data_center.proteins_dict.keys()
        for prot_key in keys_prot_dict:  # Updating the already existing protein objects functionality so that it holeds what it actually can do
            protein = self.data_center.recommendations_dict[self.data_center.proteins_dict[prot_key].name] = {}
            keys_annot_dict = self.data_center.proteins_dict[prot_key].annotations_dict.keys()
            for anno_key in keys_annot_dict:
                self.data_center.proteins_dict[prot_key].annotations_dict[anno_key].functionality.term_name = self.data_center.go_dict[anno_key]
                funct = self.data_center.proteins_dict[prot_key].annotations_dict[anno_key].functionality.term_name
                protein[funct] = float(self.data_center.proteins_dict[prot_key].annotations_dict[anno_key].evidence_code.numeric_value)


class FunctionalityInGO:  # Class for the GO objects that represent functionality of a protein
    def __init__(self):
        self.id = 0
        self.term_name = ""


class EvidenceCode:  # Class for the Evidence Code Objects that hold the acronym and the reliability score
    def __init__(self):
        self.acronym = ""
        self.numeric_value = 0


class Annotation:  # Class for the annotations of the protein
    def __init__(self):
        self.functionality = ""
        self.evidence_code = 0


class Protein:  # Class for each protein that holds its id, name and its dictionary of annotations
    def __init__(self):
        self.id = ""
        self.name = ""
        self.annotations_dict = {}


class DataCenter:  # Class that represents the data center, all of the information is hold here
    def __init__(self):
        self.proteins_dict = {}
        self.evidence_codes_dict = {}
        self.go_dict = {}
        self.recommendations_dict = {}

    def RnP_Annotation(self):  # Method for selecting GEO_human file from file browser and parsing it, while parsing creates protein objects
        self.Annotation_file = tkFileDialog.askopenfilename(initialdir="C:\Users\Yunus Stahlschmidt\PycharmProjects\ENGR 102\MP 2", title="Select Annotation File",
                                                            filetypes=(("gaf Files", "*.gaf"), ("all files", "*.*")))
        with open(self.Annotation_file, "r") as gaf_file:
            for line in gaf_file:
                if line.startswith("!"):  # skip unnecessary lines at the beginning of the file
                    continue
                else:
                    my_list = line.split('\t')
                    my_2_list = my_list[1:3]
                    my_2_list.append(my_list[4])
                    my_2_list.append(my_list[6])
                    if my_2_list[0] not in self.proteins_dict:  # checks if protein already exists in the database
                        protein = Protein()  # Creates protein objects and assigns its values from the line
                        protein.id = my_2_list[0]
                        protein.name = my_2_list[1]
                        self.proteins_dict[protein.id] = protein
                        annotation = Annotation()
                        annotation.functionality = FunctionalityInGO()
                        annotation.functionality.id = my_2_list[2]
                        annotation.evidence_code = EvidenceCode()
                        annotation.evidence_code.acronym = my_2_list[3]
                        self.proteins_dict[my_2_list[0]].annotations_dict[my_2_list[2]] = annotation
                    else:
                        annotation = Annotation()  # If the protein already exists it adds annotations to the annotation dict of the protein
                        annotation.functionality = FunctionalityInGO()
                        annotation.functionality.id = my_2_list[2]
                        annotation.evidence_code = EvidenceCode()
                        annotation.evidence_code.acronym = my_2_list[3]
                        self.proteins_dict[my_2_list[0]].annotations_dict[my_2_list[2]] = annotation

    def RnP_ECV(self):  # Method for selecting ecv file and parsing it
        self.ECV_file = tkFileDialog.askopenfilename(initialdir="C:\Users\Yunus Stahlschmidt\PycharmProjects\ENGR 102\MP 2", title="Select ECV File",
                                                     filetypes=(("txt Files", "*.txt"), ("all files", "*.*")))
        with open(self.ECV_file) as ecv_file:
            for line in ecv_file:  # Creates Evidence Code objects and stores the to the Evidence code dict of the databsse
                my_list = line.split()
                ecv_obj = EvidenceCode()
                ecv_obj.acronym = my_list[0]
                ecv_obj.numeric_value = my_list[1]
                self.evidence_codes_dict[my_list[0]] = ecv_obj

    def file_len(self):  # Looks how long the .obo file is since we have to read 2 lines at a time and dont want to take unnecessary long time for this
        with open(self.GO_file) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def RnP_GO(self):  # Method to select GO file and parse it
        self.GO_file = tkFileDialog.askopenfilename(initialdir="C:\Users\Yunus Stahlschmidt\PycharmProjects\ENGR 102\MP 2", title="Select GO File",
                                                    filetypes=(("obo Files", "*.obo"), ("all files", "*.*")))
        with open(self.GO_file) as go_file:
            for i in range(self.file_len()/2):  # Reads 2 lines at a time from the GO file and strores the info in the go dict of the database
                my_list = []
                for l in range(2):
                    my_list.append(go_file.readline()[:-1])
                my_2_list = [my_list[0][4:], my_list[1][6:]]
                self.go_dict[my_2_list[0]] = my_2_list[1]


root = Tk()
root.title("Protein Function Prediction Tool")
root.geometry("1510x600")
app = ProteinFunctionTool(root)
root.mainloop()
