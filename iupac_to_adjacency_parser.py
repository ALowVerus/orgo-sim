from molecule_builder import *
from iupac_name_parser import *
from random import randint


def name_me(mol):

    # Arbitrarily choose a root atom. This should not effect the end result.
    root = id(mol.atoms[randint(0, len(mol.atoms) - 1)])

    # Initialize items for the DFS tree and the back edges not in the tree.
    dfs_tree = {}
    back_edges = set()

    # Create a node class to contain information about the DFS tree.
    class Node:
        def __init__(self, parent, node_id, seen_count):
            self.id = node_id
            self.parent = parent
            self.seen_count = seen_count
            self.children = set()

    # Define a recursive DFS which can be used to identify back edges and thus cycles
    def dfs(atom_id, parent=None):
        dfs_tree[atom_id] = Node(parent, atom_id, len(dfs_tree))
        for neighbor_id in obj(atom_id).connections:
            if neighbor_id not in dfs_tree:
                dfs_tree[atom_id].children.add(neighbor_id)
                dfs(neighbor_id, atom_id)
            elif neighbor_id != parent:
                back_edges.add(tuple(sorted((neighbor_id, atom_id))))
    dfs(root)

    # View the results of the DFS.
    def recursively_print_dfs(atom_id, indent=0):
        if obj(atom_id).element != "H":
            print("  " * indent, str(obj(atom_id)))
        for neighbor_id in dfs_tree[atom_id].children:
            recursively_print_dfs(neighbor_id, indent + 1)
    # recursively_print_dfs(root)

    # Find the cycle parts by traversing back up through the sides of the back edges.
    def find_cycle_parts(edge):
        a, b = edge
        parts = {a, b}
        while a != b:
            if dfs_tree[a].seen_count > dfs_tree[b].seen_count:
                parts.add(dfs_tree[a].parent)
                a = dfs_tree[a].parent
            else:
                parts.add(dfs_tree[b].parent)
                b = dfs_tree[b].parent
        return parts
    cycles = []
    for back_edge in back_edges:
        cycles.append(find_cycle_parts(back_edge))
    cyclic_atoms = {}
    for c_i, cycle in enumerate(cycles):
        for atom in cycle:
            if atom not in cyclic_atoms:
                cyclic_atoms[atom] = []
            cyclic_atoms[atom].append(c_i)

    # Now that we have the cycles identified, get non-cyclic candidates for the root.
    def find_longest_chain(atom_id, parent=None, parent_count=0, indent=0):
        # If the atom at hand is in a cycle, ignore its cycle's parts and check all substituents.
        if atom_id in cyclic_atoms:
            # Create a new structure to hold new candidates.
            new_candidates = [[]]
            # Iterate over all atoms in the cycle.
            for cyclic_atom in cycles[cyclic_atoms[atom_id][0]]:
                # Iterate over all substituents of the cyclic atom.
                for neighbor_id in obj(cyclic_atom).connections:
                    # Refuse to backtrack over the initial caller's parent or the rest of the cycle.
                    if neighbor_id != parent and neighbor_id not in cycles[cyclic_atoms[atom_id][0]]:
                        # Recurse over a given substituent.
                        result = find_longest_chain(neighbor_id, cyclic_atom, 0, indent+1)
                        # For each candidate result, if it is better than or equal to the existing ones, keep it.
                        for candidate in result["New_Candidates"]:
                            # If it is better, reset.
                            if len(candidate) > len(new_candidates[0]):
                                new_candidates = [candidate]
                            # If it is equal, append.
                            elif len(candidate) == len(new_candidates[0]):
                                new_candidates.append(candidate)
                        # If the New_Atoms chain leading up to the substituent was long enough, add it too.
                        # If it is better, reset.
                        if len(result["New_Atoms"]) > len(new_candidates[0]):
                            new_candidates = [result["New_Atoms"]]
                        # If it is equal, append.
                        elif len(result["New_Atoms"]) == len(new_candidates[0]):
                            new_candidates.append(result["New_Atoms"])
            # Since we hit a cycle, and cycles are stopping points, return the list of new candidates and an empty set.
            to_return = {"New_Candidates": new_candidates, "New_Atoms": [[]]}
        # If the current atom is a Carbon, add it to all tracked chains and return the constructs upwards.
        elif obj(atom_id).element == "C":
            # Initialize new data structures to hold options.
            new_candidates = []
            new_atoms_options = [[]]
            # Gather results from all neighbors of the called atom.
            results = [find_longest_chain(neighbor_id, atom_id, parent_count+1, indent+1)
                       for neighbor_id in obj(atom_id).connections if neighbor_id != parent]
            # Check the results.
            for result in results:
                new_candidates.extend(result["New_Candidates"])
                if len(new_atoms_options) == 0 or len(result["New_Atoms"]) > len(new_atoms_options[0]):
                    new_atoms_options = result["New_Atoms"]
                elif len(result["New_Atoms"][0]) == len(new_atoms_options[0]) and len(result["New_Atoms"][0]) > 0:
                    new_atoms_options.extend(result["New_Atoms"])
            # If there is more than one result, check all results against each other for possible candidates.
            if len(results) > 1:
                for i in range(len(results) - 1):
                    for j in range(i + 1, len(results)):
                        set_a, set_b = results[i]["New_Atoms"], results[j]["New_Atoms"]
                        for a in set_a:
                            for b in set_b:
                                b = b[::-1]
                                if len(set(a).intersection(set(b))) == 0:
                                    # If the pair is better, reset.
                                    if len(new_candidates) == 0 or len(a) + len(b) + 1 > len(new_candidates[0]):
                                        new_candidates = [a + [atom_id] + b]
                                    # If the pair is equal, append.
                                    elif len(a) + len(b) + 1 == len(new_candidates[0]) and a + [atom_id] + b != new_candidates[0]:
                                        new_candidates.append(a + [atom_id] + b)
            to_return = {"New_Candidates": new_candidates, "New_Atoms": [option + [atom_id] for option in new_atoms_options]}
        # If the current is atom in neither a carbon nor a cyclic atom, ignore it, end all tracked options, and move on.
        else:
            # Initialize new data structures to hold options.
            new_candidates = []
            # Iterate over all neighbors on the called atom.
            results = [find_longest_chain(neighbor_id, atom_id, 0, indent+1)
                       for neighbor_id in obj(atom_id).connections if neighbor_id != parent]
            # Check the results.
            for result in results:
                new_candidates.extend(result["New_Candidates"])
                if len(new_candidates) == 0 or len(result["New_Atoms"]) > len(new_candidates[0]):
                    new_candidates = [result["New_Atoms"]]
                elif len(result["New_Atoms"]) == len(new_candidates[0]):
                    new_candidates.append(result["New_Atoms"])
            to_return = {"New_Candidates": new_candidates, "New_Atoms": [[]]}
        # print(" " + indent * "  " + "This atom:", str(obj(atom_id)),
        #       "\n", indent * "  " + "New_Candidates:",
        #       [[str(obj(atom_idn)) for atom_idn in atom_list] for atom_list in to_return["New_Candidates"]],
        #       "\n", indent * "  " + "New_Atoms:",
        #       [[str(obj(atom_idn)) for atom_idn in atom_list] for atom_list in to_return["New_Atoms"]])

        return to_return

    print("Longest chain:", [[str(obj(atom_idn)) for atom_idn in atom_list]
                             for atom_list in find_longest_chain(root)["New_Candidates"]])

    return "-----"


def check_out(name):
    # Convert the name to an AST
    ast = astify(tokenize(name))
    # Print the AST to confirm its worth
    recursively_print(ast)
    # Populate the molecule
    mol = Molecule(name).recursively_populate(ast).closer()
    # Parse out the formula from both methods to confirm fecundity
    name_parser_result = stringify_atoms_dict(atomize(ast))
    mol_result = mol.formula
    # Print results
    print(name_parser_result, mol_result, name_parser_result == mol_result)
    print([str(atom) for atom in mol.atoms if atom.element != "H"])
    # Run a dfs on the molecule to determine its structural hierarchy
    for i in range(20):
        mol_name = name_me(mol)
    print(mol_name)


# check_out("1-[1-[1-[methyl]methyl]methyl]methylmethylmethane")
check_out("methylmethylpropane")
