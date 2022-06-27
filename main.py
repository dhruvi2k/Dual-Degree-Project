from utils import *
from argparse import ArgumentParser
import os
# id_type = 'PMCID', id = '4916225', policy = 'one_step'
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "outputs")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)


def DDPmain(id_types, ids, policy, output_name):

    statements = []
    for id_type, id in zip(id_types, ids):
        reach_processor = processing_input(id, id_type)
        statements.extend(reach_processor.statements)

    unique_statements = getting_unique_statements(statements)
    print(unique_statements)
    model, obj1, obj2, obj3 = assembling_the_model(unique_statements, policy)
    react_list, prod_list = creating_reactant_list(obj1)
    react_components, prod_components = creating_reactant_components(
        obj1, react_list, prod_list)
    react_separate, prod_separate = creating_separate(
        react_components, prod_components)
    final_list_new = final_list(react_separate, prod_separate)
    forward_params = forward_parameters(model)
    parameter_init = parameter_init_list(forward_params)
    react_init_parameters = react_init_params(obj1, parameter_init)
    initial_l = initial_list(model)
    init_spec = init_species(parameter_init, obj1)
    react_separate_for_init, prod_separate_for_init = separate_for_init(
        react_separate, prod_separate)
    indices = get_indices(
        init_spec, react_separate_for_init, prod_separate_for_init)
    stoich_matrix = creating_stoich_matrix(
        obj1, final_list_new, react_components, prod_components)
    y = creating_output_matrix(stoich_matrix, final_list_new)
    parameter = creating_the_parameter(react_init_parameters)
    final_dict = creating_the_final_dict(
        final_list_new, react_components, prod_components, react_init_parameters)
    states = creating_the_keys(final_dict, final_list_new)
    mapped_dict = creating_a_mapped_dict(final_list_new, states)
    replace_dict(mapped_dict, final_dict)
    model_equations = generate_equations(final_dict)
    output_string, b_final, y_final = output_as_states(y, mapped_dict, states)
    first_block = creating_the_first_block(b_final)
    second_block = creating_the_second_block(parameter, states, y_final)
    third_block = creating_the_third_block(model_equations, output_string)
    fourth_block = creating_the_fourth_block()
    fifth_block = creating_the_fifth_block(
        indices, react_separate, prod_separate, mapped_dict, states, init_spec)
    final_text = creating_the_final_model(
        first_block, second_block, third_block, fourth_block, fifth_block)
    forming_the_text_file(final_text, os.path.join(OUTPUT_DIR, output_name))


if __name__ == "__main__":
    parser = ArgumentParser("indra-model")
    parser.add_argument("--id_types", type=str, nargs="+",
                        help="types of publication ID, eg PMCID")
    parser.add_argument("--ids", type=str, nargs="+",
                        help="values of publication ID, eg 4916225")
    parser.add_argument("--policy", type=str,
                        help="type of assembly policy, eg one_step")
    parser.add_argument("--output_txt", type=str,
                        help="name of the output file that needs to be generated", default="output.txt")

    args = parser.parse_args()

    DDPmain(args.id_types, args.ids, args.policy, args.output_txt)
