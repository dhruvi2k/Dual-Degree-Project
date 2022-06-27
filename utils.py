import indra
import numpy as np
from indra.sources import reach
from indra.sources import trips
from indra.assemblers.pysb import PysbAssembler
from itertools import combinations
from numpy.linalg import matrix_rank
from ast import literal_eval as make_tuple
import string
import re
from numpy.linalg import matrix_rank
import random


def processing_input(Id, Type):
    if Type == 'PMCID':
        Reach_processor = reach.process_pmc(Id)
    if Type == 'PMID':
        Reach_processor = reach.process_pubmed_abstract(Id)
    if Type == 'raw text':
        Reach_processor = reach.process_text(Id)
    return Reach_processor


def getting_the_statements(Reach_processor):
    Statements = Reach_processor.statements
    return Statements


def getting_unique_statements(Statements):
    My_dict = {}
    Unique_statements = []
    for Statement in Statements:
        if str(Statement) not in My_dict:
            Unique_statements.append((Statement))
            My_dict[str(Statement)] = 1
    return Unique_statements


def assembling_the_model(Unique_statements, Policy):
    pa = PysbAssembler()
    pa.add_statements(Unique_statements)
    model = pa.make_model(policies=Policy)
    Obj1 = pa.model.rules
    Obj2 = pa.model.monomers
    Obj3 = pa.model.parameters
    return model, Obj1, Obj2, Obj3


def creating_reactant_list(Obj1):
    React_list = []
    Prod_list = []
    for i in range(len(Obj1)):
        Curr_rule = Obj1[i]
        Curr_dict = Curr_rule.__dict__
        Curr_react = str(Curr_dict['reactant_pattern'])
        Curr_prod = str(Curr_dict['product_pattern'])
        React_list.append(Curr_react)
        Prod_list.append(Curr_prod)
    return React_list, Prod_list


def creating_reactant_components(Obj1, React_list, Prod_list):
    # making a list of reactant side species as list of lists
    React_components = list(np.zeros((len(Obj1), 2)))
    # making a list of product side species as list of lists
    Prod_components = list(np.zeros((len(Obj1), 2)))
    for i in range(len(Obj1)):
        if '+' in React_list[i]:
            React_components[i] = React_list[i].split('+')
        if '+' in Prod_list[i]:
            Prod_components[i] = Prod_list[i].split('+')
        if '+' not in React_list[i]:
            React_components[i] = [React_list[i]]
        if '+' not in Prod_list[i]:
            Prod_components[i] = [Prod_list[i]]

    for i in range(len(React_components)):
        for j in range(len(React_components[i])):
            React_components[i][j] = React_components[i][j].strip()
    for i in range(len(Prod_components)):
        for j in range(len(Prod_components[i])):
            Prod_components[i][j] = Prod_components[i][j].strip()

    return React_components, Prod_components


def creating_separate(React_components, Prod_components):
    React_separate = []
    Prod_separate = []
    # creating a list for reactants separately
    for i in range(len(React_components)):
        Curr_reac = React_components[i]
        for j in range(len(Curr_reac)):
            React_separate.append(Curr_reac[j])

    # creating a list for products
    for i in range(len(Prod_components)):
        Curr_prod = Prod_components[i]
        for j in range(len(Curr_prod)):
            Prod_separate.append(Curr_prod[j])

    return React_separate, Prod_separate


def final_list(React_separate, Prod_separate):
    Final_list_new = list(set(React_separate + Prod_separate))
    return Final_list_new


def forward_parameters(Model):
    Forward_parameters = []
    for i in range(len(Model.__dict__['parameters'])):
        Forward_parameters.append(
            str(Model.__dict__['parameters'][i]).lstrip('Parameter'))
    return Forward_parameters


def parameter_init_list(Forward_parameters):
    Parameter_init_list = []
    for i in range(len(Forward_parameters)):
        Curr_tuple = make_tuple(Forward_parameters[i])
        Curr_list = list(Curr_tuple)
        Parameter_init_list.append(Curr_list)
    return Parameter_init_list


def react_init_params(Obj1, Parameter_init_list):
    React_init_params = [[0]*2 for i in range(len(Obj1))]
    for i in range(len(Parameter_init_list)):
        for j in range(len(Obj1)):
            if Parameter_init_list[i][0] in str(Obj1.__dict__['_elements'][j]):
                React_init_params[j] = Parameter_init_list[i]
    return React_init_params


def initial_list(Model):
    Initial_list = []
    for i in range(len(Model.__dict__['initials'])):
        Initial_list.append(
            str(Model.__dict__['initials'][i]).lstrip('Initial'))

    for j in range(len(Initial_list)):
        Initial_list[j].strip()
    return Initial_list


def init_species(Parameter_init_list, Obj1):
    Number_of_init_species = len(Parameter_init_list) - len(Obj1)
    Init_species = Parameter_init_list[-Number_of_init_species:]
    for i in range(len(Init_species)):
        Init_species[i][0] = Init_species[i][0].replace("_0", "")
    return Init_species


def separate_for_init(React_separate, Prod_separate):
    React_separate_for_init = []
    Prod_separate_for_init = []
    for i in range(len(React_separate)):
        React_separate_for_init.append(
            re.sub('\(.*?\)', "", React_separate[i]))
    for j in range(len(Prod_separate)):
        Prod_separate_for_init.append(re.sub('\(.*?\)', "", Prod_separate[j]))
    return React_separate_for_init, Prod_separate_for_init


def get_indices(Init_species, React_separate_for_init, Prod_separate_for_init):
    Indices = []
    for i in range(len(Init_species)):
        if Init_species[i][0] in React_separate_for_init:
            Num = -React_separate_for_init.index(Init_species[i][0])-1
            Indices.append(Num)

        else:
            Num = Prod_separate_for_init.index(Init_species[i][0])+1
            Indices.append(Num)

    return Indices


def creating_stoich_matrix(Obj1, Final_list_new, React_components, Prod_components):
    Stoich_matrix = np.zeros((len(Obj1), len(Final_list_new)))
    for i in range(Stoich_matrix.shape[1]):
        for j in range(Stoich_matrix.shape[0]):
            # here assuming that if the component appears on the
            if (Final_list_new[i] in React_components[j]):
                Stoich_matrix[j, i] = -1
                break
            if (Final_list_new[i] in Prod_components[j]):
                Stoich_matrix[j, i] = 1
                break
            else:
                Stoich_matrix[j, i] = 0

    return Stoich_matrix


def creating_output_matrix(Stoich_matrix, Final_list_new):
    Y = []
    N = matrix_rank(Stoich_matrix)
    if Stoich_matrix.shape[1] > 20:
        random_list = []
        for i in range(10):
            Curr_list = random.sample(range(Stoich_matrix.shape[1]), N)
            random_list.append(Curr_list)

        for i in random_list:
            if np.linalg.matrix_rank(Stoich_matrix[:, i]) == N:
                Curr_list = i
                break

    else:
        Perm_list = []
        for i in range((Stoich_matrix.shape[1])):
            Perm_list.append(i)
        for i in combinations(Perm_list, N):
            Curr_list = list(i)
            Curr_matrix = Stoich_matrix[:, Curr_list]
            if matrix_rank(Curr_matrix) == N:
                break

    for i in range(len(Curr_list)):
        Y.append(Final_list_new[i])

    return Y


def creating_the_parameter(React_init_params):
    Parameter = []
    for i in range(len(React_init_params)):
        Parameter.append(React_init_params[i][0])
    return Parameter


def creating_the_final_dict(Final_list_new, React_components, Prod_components,  React_init_params):
    Parameter = []
    for i in range(len(React_init_params)):
        Parameter.append(React_init_params[i][0])
    Final_dict = dict([(key, []) for key in Final_list_new])
    for i in range(len(Final_list_new)):
        to_append_final = []
        for j in range(len(React_components)):
            if Final_list_new[i] in React_components[j]:
                to_append_list = [-1, Parameter[j]]
                React_list = React_components[j]
                for k in range(len(React_list)):
                    to_append_list.append(React_list[k])
                to_append_final.append(to_append_list)
            if (Final_list_new[i] in Prod_components[j]) and (Final_list_new[i] not in React_components[j]):
                to_append_list = [1, Parameter[j]]
                React_list = React_components[j]
                for k in range(len(React_list)):
                    to_append_list.append(React_list[k])
                to_append_final.append(to_append_list)

    # to_append_list.append(to_append_list)
        Final_dict[Final_list_new[i]] = to_append_final
    return Final_dict


def creating_the_keys(Final_dict, Final_list_new):
    States = []
    for i in range(len(Final_list_new)):
        Curr = 'x' + str(i+1)
        States.append(Curr)
    return States


def creating_a_mapped_dict(Final_list_new, States):
    Mapped_dict = {Final_list_new[i]: States[i] for i in range(len(States))}
    return Mapped_dict


def replace_dict(d1, d2):
    for k, v in d2.items():

        for ls in v:
            for i in range(len(ls)):
                if ls[i] in d1:
                    ls[i] = d1[ls[i]]

    for k in d1:
        if k in d2:
            d2[d1[k]] = d2.pop(k)


def applying(Mapped_dict, Final_dict):
    Mapped_dict = replace_dict(Mapped_dict, Final_dict)
    return Mapped_dict


def generate_equations(Final_dict):
    Model_equations = ""
    for k, v in Final_dict.items():
        eqn = f"df({k}, t)=" + \
            " + ".join(["*".join(list(map(str, ls))) for ls in v])
        Model_equations += eqn + ","
    Model_equations = Model_equations.strip(",")
    return Model_equations


def output_as_states(Y, Mapped_dict, States):
    Output_string = ""
    Y_final = []
    Y_as_states = []
    for i in range(len(Y)):
        Y_final.append('y' + str(i+1))

    for i in range(len(Y)):
        Y_as_states.append(Mapped_dict[Y[i]])

    for i in range(len(Y)):
        Output_string = Output_string + \
            ",  %s = %s" % (Y_final[i], Y_as_states[i])

    B_final = Y_final + States
    return Output_string, B_final, Y_final


def creating_the_first_block(B_final):
    New_string = str(B_final).replace("'", "")
    New_string1 = New_string.replace("[", "{")
    New_string2 = New_string1.replace("]", "}")
    B = 'B_:=' + New_string2 + '$'
    First_block = ('''WRITE "FINAL_MODEL"$''' + "\n" + B +
                   "\n" + "FOR EACH EL_ IN B_ DO DEPEND EL_, T$")
    return First_block


def creating_the_second_block(Parameter, States, Y_final):
    Parameter = str(Parameter).replace("'", "")
    Parameter = Parameter.replace("[", "{")
    Parameter = Parameter.replace("]", "}")
    Parameter_string = "B1_:=" + Parameter + "$"
    Second_block = (Parameter_string + "\n" + "NX_:=%s$" % len(States) +
                    "\n" + "NU_:= 0$" + "\n" + "NY_:= %s$" % len(Y_final))
    return Second_block


def creating_the_third_block(Model_equations, Output_string):
    Third_block = "C_:= {" + Model_equations + Output_string + "}$"
    return Third_block


def creating_the_fourth_block():
    Fourth_block = ("FLAG_:=1$" "\n" + "B2_:={}$")
    return Fourth_block


def creating_the_fifth_block(Indices, React_separate, Prod_separate, Mapped_dict, States, Init_species):
    States_to_be_included = []
    for i in range(len(Indices)):
        index = Indices[i]
        if index < 0:
            index = -1*index - 1
            dum = React_separate[index]
            States_to_be_included.append(Mapped_dict[dum])
        else:
            index = index - 1
            dum = Prod_separate[index]
            States_to_be_included.append(Mapped_dict[dum])

    States_not_included = list(set(States) - set(States_to_be_included))
    main1 = "%s" % States_not_included[0]
    for i in range(1, len(States_not_included)-1):
        main1 = main1 + "," + "%s" % States_not_included[i]
    n = len(States_not_included)-1
    main1 = "{" + main1 + ",%s}" % States_not_included[n] + "$"

    main = "%s = %s" % (States_to_be_included[0], Init_species[0][1])
    for i in range(1, len(States_to_be_included)):
        main = main + ", " + \
            "%s = %s" % (States_to_be_included[i], Init_species[i][1])

    Fifth_block = ("daisy()$" + "\n" + "ICK_:= {" + main + "}$" +
                   "\n" + "ICUNK_:=%s" % main1 + "\n" + "CONDINIZ$" + "\n" + "END$")

    return Fifth_block


def creating_the_final_model(First_block, Second_block, Third_block, Fourth_block, Fifth_block):
    Final_text = (First_block + "\n" + Second_block + "\n" +
                  Third_block + "\n" + Fourth_block + "\n" + Fifth_block)
    return Final_text


def forming_the_text_file(Final_text, output_path):
    with open(output_path, "w") as text_file:
        text_file.write("%s" % Final_text)
