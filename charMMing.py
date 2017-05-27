__author__ = "kbuzar"

from os import makedirs, path, walk
from subprocess import *
import re
try:
    import readline
except ImportError:
    import pyreadline as readline

class MyCompleter(object):  # Custom completer

    def __init__(self, options):
        self.options = sorted(options)

    def complete(self, text, state):
        if state == 0:  # on first trigger, build possible matches
            if text:  # cache matches (entries that start with entered text)
                self.matches = [s for s in self.options if s and s.startswith(text)]
            else:  # no text entered, all matches possible
                self.matches = self.options[:]

        # return match indexed by state
        try:
            return self.matches[state]
        except IndexError:
            return None


class CharMMing:

    def __init__(self, pdb_file):
        #self.debug_mode = True
        self.root_dir = "/home/kbuzar/charmm-files/data"
        self.charmm_path = ""
        self.pdb_file = pdb_file
        self.create_directory(self.pdb_file[:-4])
        self.pdb_dir = "./%s/pdb" % self.pdb_file[:-4]
        self.create_directory(self.pdb_dir)
        self.coord_dir = "./%s/coord" % self.pdb_file[:-4]
        self.create_directory(self.coord_dir)
        self.aminoacid_lenght_dict = {"ALA": 5, "LEU": 8, "ASP": 8, "PRO": 7, "PHE": 11,
                                      "ASN": 8, "HIS": 10, "HSE": 10, "ARG": 11, "GLU": 9,
                                      "MET": 8, "VAL": 7, "ILE": 8, "TRP": 14, "TYR": 12,
                                      "GLY": 4, "SEC": 6, "CYS": 6, "THR": 7, "SER": 6,
                                      "GLN": 9, "LYS": 9
                                      }

    def create_directory(self, directory_path):
        if not path.exists(directory_path):
            makedirs(directory_path)

    def read_PDB(self):
        with open(self.pdb_file, "r") as file:
            file_lines = [line.strip() for line in file]
            self.seq_dict = {}
            self.seq_aminoacid_number_dict = {}
            self.het_dict = {}
            self.coord_dict = {}
            self.missing_dict = {}
            self.missing_dict_values = {}
            for theline in file_lines:
                theline = theline.split()
                if theline[0] == "SEQRES":
                    curr_letter = theline[2]
                    self.seq_dict[curr_letter] = []
                    self.seq_aminoacid_number_dict[theline[2]] = int(theline[3])
                    self.coord_dict[curr_letter] = []
                    self.missing_dict[curr_letter] = []
                    self.missing_dict_values[curr_letter] = 0
                # elif theline[0] == "HET":
                #     self.het_dict[theline[1]] = []
                #     self.coord_dict[theline[1]] = []
            for aline in file_lines:
                aline = aline.split()
                if aline[0] == "SEQRES":
                    seq_part = aline[4:]
                    for aminoacid in seq_part:
                        self.seq_dict[aline[2]].append(aminoacid)
                #elif aline[0] == "HET":
                #    pass
                elif aline[0] == "ATOM":
                    try:
                        curr_letter = aline[4]
                        self.coord_dict[curr_letter].append(aline)
                    except KeyError:
                        curr_letter = aline[3]
                        self.coord_dict[curr_letter].append(aline)
                elif aline[0] == "TER":
                    curr_letter = aline[3]
                #     self.coord_dict[curr_letter].append(aline)
                # elif aline[0] == "HETATOM":
                #     pass
            for bline in file_lines:
                bline = bline.split()
                if bline[0] == "REMARK" and bline[1] == "465" and len(bline) == 5:
                    self.missing_dict[bline[3]].append(int(bline[4])-1)
                    self.missing_dict_values[bline[3]] += 1
        for key, value in self.seq_dict.iteritems():
            increment_value = 0
            for element in self.missing_dict[key]:
                try:
                    del self.seq_dict[key][element-increment_value]
                except IndexError:
                    del self.seq_dict[key][-1]
                increment_value += 1
            self.seq_aminoacid_number_dict[key] -= self.missing_dict_values[key]
        self.cut_pdb()

    def cut_pdb(self):
        for key, value in self.coord_dict.iteritems():
            full_path = "./%s/pdb/chain_%s.pdb" % (self.pdb_file[:-4], key)
            print full_path
            with open(full_path, "w") as f:
                aminoacid_counter = 1
                emergency_print_list = []
                general_counter = 0
                line_counter = 0
                lenght_of = len(value)
                for line in value:
                    line_counter += 1
                    if len(line) != 12:
                        if len(line[2]) > 3:
                            amino = line[2][-4:]
                            atom = line[2][:-4]
                            line = line[:2] + [atom] + [amino] + line[2:]
                            del line[4]
                        if line[-2] > 5:
                            if len(line[-2]) > 5:
                                value = line[-2][:4]
                                rest = line[-2][4:]
                                line = line[:-2] + [value] + [rest] + line[-2:]
                                del line[-2]
                        if len(line[6]) > 7:
                            test_case = line[6].split("-")
                            if len(test_case) == 3 and len(test_case[0]) > 0:
                                first = test_case[0]
                                second = "-" + test_case[1]
                                third = "-" + test_case[2]
                                line = line[:5] + [first] + [second] + [third] + line[5:]
                                del line[8]
                            elif len(test_case) == 3:
                                first = "-" + test_case[1]
                                second = "-" + test_case[2]
                                line = line[:5] + [first] + [second] + line[5:]
                                del line[7]
                            elif len(test_case) == 4:
                                first = "-" + test_case[0]
                                second = "-" + test_case[1]
                                third = "-" + test_case[2]
                                line = line[:5] + [first] + [second] + [third] + line[5:]
                                del line[8]
                            elif len(test_case) == 2:
                                first = test_case[1]
                                second = "-" + test_case[2]
                                line = line[:5] + [first] + [second] + line[5:]
                                del line[7]
                        elif len(line[7]) > 7:
                            test_case = line[6].split("-")
                            if len(test_case) == 3:
                                first = "-" + test_case[1]
                                second = "-" + test_case[2]
                                line = line[:6] + [first] + [second] + line[6:]
                            elif len(test_case) == 2:
                                first = test_case[0]
                                second = "-" + test_case[1]
                                line = line[:6] + [first] + [second] + line[6:]
                                del line[8]
                    if len(line[3]) == 4:
                        if line[3][0] == "B":
                            general_counter -= 1
                            continue
                        else:
                            line[3] = line[3][1:]
                    if line[3] == "HIS":
                        line[3] = "HSE"
                    if len(emergency_print_list) > 0 and general_counter != 0:
                        if emergency_print_list[-1][3] != line[3]:
                            limiter = self.aminoacid_lenght_dict[emergency_print_list[-1][3]]
                            counter = 0
                            for entriees in emergency_print_list:
                                entriees = self.line_creator(entriees, key, aminoacid_counter)
                                print >> f, entriees
                            # if len(emergency_print_list) == limiter or len(emergency_print_list) % limiter == 0:
                            #     for entries in emergency_print_list:
                            #         if counter == limiter:
                            #             counter = 0
                            #             aminoacid_counter += 1
                            #         else:
                            #             entries = self.line_creator(entries, key, aminoacid_counter)
                            #             print >> f, entries
                            #             counter += 1
                            # else:
                            #     for entriees in emergency_print_list:
                            #         entriees = self.line_creator(entriees, key, aminoacid_counter)
                            #         print >> f, entriees
                            #         aminoacid_counter += 1
                            emergency_print_list = []
                        if lenght_of == line_counter:
                            emergency_print_list.append(line)
                            for entries in emergency_print_list:
                                entries = self.line_creator(entries, key, aminoacid_counter)
                                print >> f, entries
                    emergency_print_list.append(line)
                    general_counter += 1
                print >>f, "END\n\n"

    def line_creator(self, line, key, aminoacid_counter):
        prox_identifier = "PRO%s" % key
        print_line = line[0]
        print_line += " " * (7-len(line[1]))
        print_line += line[1]
        if len(line[2]) == 1:
            print_line += "  %s   %s %s" % (line[2], line[3], line[4])
        elif len(line[2]) == 2:
            print_line += "  %s  %s %s" % (line[2], line[3], line[4])
        elif len(line[2]) == 3 and len(line[3]) != 4:
            print_line += " %s  %s %s" % (line[2], line[3], line[4])
        else:
            print_line += " %s %s %s" % (line[2], line[3], line[4])
        print_line += " " * (4-len(str(aminoacid_counter)))
        print_line += str(aminoacid_counter)
        print_line += " " * (12-len(line[6]))
        print_line += line[6]
        print_line += " " * (8 - len(line[7]))
        print_line += line[7]
        print_line += " " * (8 - len(line[8]))
        print_line += line[8]
        print_line += "  %s %s" % (line[9], line[10])
        print_line += "      %s %s" % (prox_identifier, line[11])
        return print_line

    def create_charmm_inp_1(self):
        path = "./%s/" % self.pdb_file[:-4]
        for key, value in self.seq_aminoacid_number_dict.iteritems():
            cor_1_path = "chain_%s_read_1.cor" % key
            pdb_path = "chain_%s.pdb" % key
            with open(path + "read_chain_%s_1.inp" % key, "w") as file:
                print >>file, "* read chain A from %s" % self.pdb_file[:-4]
                print >>file, "*\n"

                print >>file, "prnlev 3"
                print >>file, "bomlev -1\n"

                print >>file, "set ROOT %s" % self.root_dir
                print >>file, "open read card unit 1 name @ROOT/top_all36_prot.rtf"
                print >>file, "read RTF card unit 1"
                print >>file, "close unit 1\n"

                print >>file, "open read card unit 1 name @ROOT/par_all36_prot.prm"
                print >>file, "read PARAM card unit 1"
                print >>file, "close unit 1\n"

                print >>file, "read sequence card"
                print >>file, "* sequence file"
                print >>file, "*"
                print >>file, "%s" % value

                aminoacid_counter = 1
                for aminoacid in self.seq_dict[key]:
                    if aminoacid == "HIS":
                        aminoacid = "HSE"
                    if aminoacid_counter <= 12:
                        print >> file, aminoacid,
                    else:
                        print >>file, aminoacid
                        aminoacid_counter = 1
                    aminoacid_counter += 1

                print >>file, "\ngenerate proa first nter last cter setup\n"

                print >>file, "open read card unit 17 name pdb/%s" % pdb_path
                print >>file, "read coor pdb unit 17"
                print >>file, "close unit 17\n"

                print >>file, "delete atom sele hydrogen end\n"

                print >>file, "ic para"
                print >>file, "ic build\n"

                print >>file, "rename resn HSD sele ((resi 205 .or. resi 262) .and. resn HSE) end\n"

                print >>file, "open write card unit 17 name coord/%s" % cor_1_path
                print >>file, "write coor card unit 17\n"

                print >>file, "stop\n"


    def create_charmm_inp_2(self):
        path = "./%s/" % self.pdb_file[:-4]
        for key, value in self.seq_aminoacid_number_dict.iteritems():
            cor_1_path = "chain_%s_read_1.cor" % key
            cor_2_path = "chain_%s_read_2.cor" % key
            with open(path + "read_chain_%s_2.inp" % key, "w") as file:
                print >>file, "* read chain A from %s" % self.pdb_file[:-4]
                print >>file, "*\n"

                print >>file, "prnlev 3"
                print >>file, "bomlev -1\n"

                print >>file, "set ROOT %s" % self.root_dir
                print >>file, "open read card unit 1 name @ROOT/top_all36_prot.rtf"
                print >>file, "read RTF card unit 1"
                print >>file, "close unit 1\n"

                print >>file, "open read card unit 1 name @ROOT/par_all36_prot.prm"
                print >>file, "read PARAM card unit 1"
                print >>file, "close unit 1\n"

                print >>file, "read sequence card"
                print >>file, "* sequence file"
                print >>file, "*"
                print >>file, "%s" % value

                aminoacid_counter = 1
                for aminoacid in self.seq_dict[key]:
                    if aminoacid == "HIS":
                        aminoacid = "HSE"
                    if aminoacid_counter <= 12:
                        print >> file, aminoacid,
                    else:
                        print >>file, aminoacid
                        aminoacid_counter = 1
                    aminoacid_counter += 1

                print >>file, "\ngenerate proa first nter last cter setup\n"

                print >>file, "open read card unit 17 name coord/%s" % cor_1_path
                print >>file, "read coor pdb unit 17"
                print >>file, "close unit 17\n"

                print >>file, "hbuild\n"

                print >>file, "open write card unit 17 name coord/%s" % cor_2_path
                print >>file, "write coor card unit 17\n"

                print >>file, "stop\n"

    def execution_interator(self):
        pass

    def execute_charmm_inp_1(self, input_path):
        pass
        execute = Popen("%s < %s" % (self.charmm_path, input_path), stderr=PIPE, stdout=PIPE, stdin=PIPE,
                        shell=True)#, cwd=path.join(root))
        err, out = execute.communicate()


    def execute_charmm_inp_2(self, input_path):
        pass
        execute = Popen("%s < %s" % (self.charmm_path, input_path), stderr=PIPE, stdout=PIPE, stdin=PIPE,
                        shell=True)#, cwd=path.join(root))
        err, out = execute.communicate()


if __name__ == '__main__':
    print(chr(27) + "[2J")  # clear screen
    print ">>>>>>>>>>>>>>>>> PRINCE CHARMMING v0.0.1<<<<<<<<<<<<<<<<\n"
    menu = True
    debug_mode = True
    filename_list = []
    filter = ".*.pdb"
    pattern = re.compile(filter)
    base_paths = walk('./').next()[2]
    for files in base_paths:
        if pattern.match(files):
            filename_list.append(files)
    while menu != False:
        print(chr(27) + "[2J")  # clear screen
        print "My Lord, should You: "
        print "[1] Update and generate STATUS report."
        print "[0] Quit thy quest and doom us all."
        try:
            #x = raw_input("What is thy bidding Master: ")
            x=1
            x = int(x)
        except SyntaxError:
            print(chr(27) + "[2J")  # clear screen
            print "\nSire, no value was entered. I beg You, reconsider thou request.\n"
            continue
        except NameError:
            print(chr(27) + "[2J")  # clear screen
            print "\nSire, Thou shalt enter only teh integers.\n"
            continue
        except ValueError:
            print(chr(27) + "[2J")  # clear screen
            print "\nSire, Thou shalt enter only teh integers.\n"
            continue
        if x == 1:
            print(chr(27) + "[2J")  # clear screen
            print "\nAvaiable PDB structures:"
            for file in filename_list:
                print file,
            print "\n"
            if debug_mode:
                pdb_file = "3wu2.pdb"
            else:
                completer = MyCompleter(filename_list)
                readline.set_completer(completer.complete)
                readline.parse_and_bind('tab: complete')
                pdb_file = raw_input("Specify the PDB file to import and CHARMM: ")
            C = CharMMing(pdb_file)
            C.read_PDB()
            C.create_charmm_inp_1()
            C.create_charmm_inp_2()
            menu = False
        elif x == 0:
            print(chr(27) + "[2J")  # clear screen
            print "End-time comes."
            print "Everybody dies."
            menu = False
        else:
            print(chr(27) + "[2J")  # clear screen
            print "That is not right sire. I beg You Master, thou shalt properly specify the position in the menu. \n"
            continue