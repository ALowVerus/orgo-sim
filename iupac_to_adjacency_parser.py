from iupac_name_parser import *
from random import randint


# Create a node class to contain information about the DFS tree.
class DFSNode:
    def __init__(self, parent, node_id, seen_count):
        self.id = node_id
        self.parent = parent
        self.seen_count = seen_count
        self.children = set()


# # Create a node class and an edge to contain information about a reduced tree
# class ReducedNode:
#     def __init__(self, classification, parts, edges):
#         self.type = classification
#         self.parts = parts
#         self.edges = edges
#
#
# class ReducedEdge:
#     def __init__(self, front, end, parts):
#         self.front = front
#         self.end = end
#         self.parts = parts


# Run a DFS to find cycles and organize the molecule as a reduced graph
def dfs_procedure(root):
    # Initialize items for the DFS tree and the back edges not in the tree.
    dfs_tree = {}
    back_edges = set()

    # Define a recursive DFS which can be used to identify back edges and thus cycles
    def dfs(atom_id, parent=None):
        dfs_tree[atom_id] = DFSNode(parent, atom_id, len(dfs_tree))
        for neighbor_id in obj(atom_id).connections:
            if neighbor_id not in dfs_tree:
                dfs_tree[atom_id].children.add(neighbor_id)
                dfs(neighbor_id, atom_id)
            elif neighbor_id != parent:
                back_edges.add(tuple(sorted((neighbor_id, atom_id))))

    dfs(root)

    # Find the cycle parts by traversing back up through the sides of the back edges.
    cycles = []
    for back_edge in back_edges:
        a, b = back_edge
        a_visited = []
        b_visited = []
        while a != b:
            if dfs_tree[a].seen_count > dfs_tree[b].seen_count:
                a_visited.append(a)
                a = dfs_tree[a].parent
            else:
                b_visited.append(b)
                b = dfs_tree[b].parent
        total_visited = a_visited + [a] + b_visited[::-1]
        cycles.append(total_visited)

    # Add the index of each cycle to a list associated with each atom,
    # so that we can list all the cycles that an atom is in.
    cyclic_atoms = {}
    for c_i, cycle in enumerate(cycles):
        for atom in cycle:
            if atom not in cyclic_atoms:
                cyclic_atoms[atom] = []
            cyclic_atoms[atom].append(c_i)

    # # Generate a reduced graph, using cycles and multi-neighbor nodes as joints and lists as edges
    # chunks = []
    # seen_atoms = set()
    # for atom_id in cyclic_atoms:
    #     if atom_id not in seen_atoms:
    #         next_chunk = {"atoms": set(), "cycles": set()}
    #         q = set(cyclic_atoms[atom_id])
    #         while q:
    #             cycle_number = q.pop()
    #             next_chunk["cycles"].add(cycle_number)
    #             for cyclic_atom_id in cycles[cycle_number]:
    #                 if cyclic_atom_id not in seen_atoms:
    #                     seen_atoms.add(cyclic_atom_id)
    #                     next_chunk["atoms"].add(cyclic_atom_id)
    #                     for cycle_id in cyclic_atoms[cyclic_atom_id]:
    #                         if cycle_id not in next_chunk:
    #                             q.add(cycle_id)
    #         chunks.append(next_chunk)
    # atoms_to_chunks = {}
    # for i, chunk in chunks:
    #     for atom_id in chunk["atoms"]:
    #         atoms_to_chunks[atom_id] = i
    #
    # # Generate a reduced graph
    # def recursively_reduce(curr_node, parent_path):
    #     if curr_node in atoms_to_chunks:
    #         chunk_atoms = chunks[atoms_to_chunks[curr_node.id]]["atoms"]
    #         reduced_node = ReducedNode("chunk", chunk_atoms, [])
    #         if parent_path:
    #             parent_path.end = reduced_node
    #             reduced_node.
    #
    #
    # reduced_root = recursively_reduce(root, None)

    # Return self-contained product objects
    return cycles, cyclic_atoms


# Given a molecule object, get its name
def find_chains(mol):

    # Arbitrarily choose a root atom. This should not effect the end result.
    root = id(mol.atoms[randint(0, len(mol.atoms) - 1)])

    cycles, cyclic_atoms = dfs_procedure(root)

    # Now that we have the cycles identified, get non-cyclic candidates for the root.
    def find_longest_chain(atom_id, parent=None, parent_count=0, indent=0):
        # Initialize new data structures to hold options.
        new_candidates = set()
        new_atoms_options = set()
        # If the atom at hand is in a cycle, ignore its cycle's parts and check all substituents.
        if atom_id in cyclic_atoms:
            # Set the cycle as a new candidate
            new_candidates.add(tuple(cycles[cyclic_atoms[atom_id][0]]))
            # Iterate over all atoms in the cycle to get their substituents.
            for cyclic_atom in cycles[cyclic_atoms[atom_id][0]]:
                # Iterate over all substituents of the cyclic atom.
                for neighbor_id in obj(cyclic_atom).connections:
                    # Refuse to backtrack over the initial caller's parent or the rest of the cycle.
                    if neighbor_id != parent and neighbor_id not in cycles[cyclic_atoms[atom_id][0]]:
                        # Recurse over a given substituent.
                        sub_chain_result = find_longest_chain(neighbor_id, cyclic_atom, 0, indent+1)
                        # For each candidate result, if it is better than or equal to the existing ones, keep it.
                        for candidate in sub_chain_result["New_Candidates"]:
                            # If it is better, reset.
                            if len(candidate) > len(new_candidates[0]):
                                new_candidates = [candidate]
                            # If it is equal, append.
                            elif len(candidate) == len(new_candidates[0]):
                                new_candidates.append(candidate)
                        # If the New_Atoms chain leading up to the substituent was long enough, add it too.
                        # If it is better, reset.
                        if len(new_candidates) == 0 \
                                or len(sub_chain_result["Unaffiliated_Chains"]) > len(next(iter(new_candidates))):
                            new_candidates = sub_chain_result["Unaffiliated_Chains"]
                        # If it is equal, append.
                        elif len(sub_chain_result["Unaffiliated_Chains"]) == len(next(iter(new_candidates))):
                            new_candidates |= sub_chain_result["Unaffiliated_Chains"]
        # If the current atom is a Carbon, add it to all tracked chains and return the constructs upwards.
        elif obj(atom_id).element == "C":
            # Gather results from all neighbors of the called atom.
            results = [find_longest_chain(neighbor_id, atom_id, parent_count+1, indent+1)
                       for neighbor_id in obj(atom_id).connections if neighbor_id != parent]
            # Check the results.
            for result in results:
                new_candidates |= result["New_Candidates"]
                if len(new_atoms_options) == 0 or len(result["Unaffiliated_Chains"]) > len(next(iter(new_atoms_options))):
                    new_atoms_options = result["Unaffiliated_Chains"]
                else:
                    new_atoms_options |= result["Unaffiliated_Chains"]
            # If there is more than one result, check all results against each other for possible candidates.
            if len(results) > 1:
                for i in range(len(results) - 1):
                    for j in range(i + 1, len(results)):
                        set_a, set_b = results[i]["Unaffiliated_Chains"], results[j]["Unaffiliated_Chains"]
                        for a in set_a:
                            for b in set_b:
                                # print(new_candidates, len(new_candidates), len(a), len(b))
                                if len(set(a).intersection(set(b))) == 0:
                                    b = b[::-1]
                                    # If the pair is better, reset.
                                    if len(new_candidates) == 0 or len(a) + len(b) + 1 > len(next(iter(new_candidates))):
                                        new_candidates = {tuple(a + (atom_id,) + b)}
                                    # If the pair is equal, append.
                                    elif len(a) + len(b) + 1 == len(next(iter(new_candidates))) \
                                            and tuple(set(a) | {atom_id} | set(b)) != set(next(iter(new_candidates))):
                                        new_candidates.add(tuple(a + (atom_id,) + b))
            # Plop this carbon onto the working atom lists if they exist
            for option in new_atoms_options.copy():
                new_atoms_options.remove(option)
                new_atoms_options.add(option + (atom_id,))
            # Generate a new atom list if this is the far end
            if len(new_atoms_options) == 0:
                new_atoms_options.add((atom_id,))
        # If the current is atom in neither a carbon nor a cyclic atom, ignore it, end all tracked options, and move on.
        else:
            # Iterate over all neighbors on the called atom.
            results = [find_longest_chain(neighbor_id, atom_id, 0, indent+1)
                       for neighbor_id in obj(atom_id).connections if neighbor_id != parent]
            # Check the results.
            for result in results:
                print(new_candidates, result["New_Candidates"])
                new_candidates |= result["New_Candidates"]
                if len(new_candidates) == 0 or len(result["Unaffiliated_Chains"]) > len(next(iter(new_candidates))):
                    new_candidates = result["Unaffiliated_Chains"]
                elif len(next(iter(new_candidates))) == len(next(iter(new_candidates))):
                    new_candidates |= result["Unaffiliated_Chains"]
        return {"New_Candidates": new_candidates, "Unaffiliated_Chains": new_atoms_options}

    def reconcile_last_chain(root):
        # Get options
        result = find_longest_chain(root)
        # Consider the new candidates as well
        print(result["New_Candidates"], result["Unaffiliated_Chains"])
        t_result = result["New_Candidates"] | result["Unaffiliated_Chains"]
        # Remove candidates with too few nodes
        t_result = [part for part in t_result if len(part) == max([len(sub_l) for sub_l in t_result])]
        # Remove cyclical duplicates
        sifted_res = []
        seen = set()
        for part in t_result:
            candy = tuple(sorted(part))
            if candy not in seen:
                sifted_res.append(part)
                seen.add(candy)
        t_result = sifted_res
        # Print result
        print("RES:", len(t_result[0]), t_result)
        for i in range(len(t_result[0])):
            t = t_result[0][i]
            print("\t", obj(t))

        # Add the New Atoms
        result = [[str(obj(atom_id)) for atom_id in tup] for tup in result["New_Candidates"] | result["Unaffiliated_Chains"]]
        # Parse out things that are too small
        result = [part for part in result if len(part) == max([len(sub_l) for sub_l in result])]
        # Parse out repeats
        result = set([tuple(sorted(part)) for part in result])
        return result

    long_chains = list(reconcile_last_chain(root))
    return long_chains


def astify_chains(long_chains):
    working_chains = long_chains + long_chains[::-1]


def resolve_true_ast(chain_asts):
    return None


def name_mol(mol):
    # Get the longest root chains
    long_chains = find_chains(mol)
    # Astify the root chains
    chain_asts = astify_chains(long_chains)
    # Determine a winning ast, if needed
    true_ast = resolve_true_ast(chain_asts)
    # Resolve into a name
    mol_name = resolve_true_ast(true_ast)
    # Return your result!
    return mol_name


def check_out(name):
    # Convert the name to an AST
    ast = astify(tokenize(name))
    # Print the AST to confirm its worth
    recursively_print(ast)
    # Populate the molecule
    mol = Molecule(name).recursively_populate(ast).closer()
    # Parse out the formula from both methods to confirm fecundity
    name_parser_result = stringify_atoms_dict(ast.atomize())
    mol_result = mol.formula
    # Print results
    print(name_parser_result, mol_result, name_parser_result == mol_result)
    # Name the molecule algorithmically
    print("long_chains:", name_mol(mol))



astify(tokenize("methylmethandiol"))
astify(tokenize("2,2-dimethylpropane"))

check_out("1-[1-[1-[methyl]methyl]methyl]methylmethylmethane")
check_out("2,2-dimethylpropane")
check_out("butane")
check_out("benzene")
check_out("2-methylbutane")
check_out("4-methyl-3-butyl-2-methylundecane")
