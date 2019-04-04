from collections import deque, defaultdict

from glypy.algorithms.similarity import commutative_similarity
from glypy.utils import root


class TopologicalInclusionMatcher(object):
    def __init__(self, target, reference, substituents=True, tolerance=0, visited=None):
        self.target = target
        self.reference = reference
        self.substituents = substituents
        self.tolerance = tolerance
        self.visited = visited or set()

    @classmethod
    def compare(cls, target, reference, substituents=True, tolerance=0, visited=None):
        inst = cls(
            target, reference, substituents=substituents,
            tolerance=tolerance, visited=visited)
        score = inst.test()
        return score

    def test_similarity(self):
        similar = commutative_similarity(
            self.target, self.reference, self.tolerance,
            include_substituents=self.substituents)
        return similar

    def test_children(self):
        match_index = self._build_child_pair_score_map()
        required_nodes = {node.id for p, node in self.target.children()}
        optimal_assignment, score = self.optimal_assignment(match_index, required_nodes)
        return optimal_assignment, score, bool(required_nodes)

    def test(self):
        key = (self.target.id, self.reference.id)
        if key in self.visited:
            return True
        self.visited.add(key)
        similar = self.test_similarity()
        if similar:
            child_pairs, score, had_children = self.test_children()
            if child_pairs:
                return similar + score
            elif not had_children:
                return similar
            else:
                return 0
        return 0

    def _build_child_pair_score_map(self):
        combinations = dict()
        for a_pos, a_child in self.target.children():
            for b_pos, b_child in self.reference.children():
                inclusion_score = self.compare(
                    a_child, b_child, substituents=self.substituents,
                    tolerance=self.tolerance, visited=self.visited)
                if inclusion_score > 0:
                    combinations[a_child.id, b_child.id] = inclusion_score
        return combinations

    def build_unique_index_pairs(self, pairs):
        '''
        Generate all unique non-overlapping sets of pairs, given in
        `pairs`
        '''
        depth = 0
        pairings = defaultdict(set)
        for a, b in pairs:
            pairings[a].add(b)
        next_current = [()]
        options_a = list(pairings)
        partial_solutions = set()
        while depth < len(options_a):
            for current in next_current:
                if len(current) > 0:
                    current_a, current_b = map(set, zip(*current))
                else:
                    current_b = set()
                components = set()
                a = options_a[depth]
                for b in pairings[a] - current_b:
                    components.add((current + ((a, b),)))
                partial_solutions.update(components)
            depth += 1
            next_current = partial_solutions
            partial_solutions = set()
        return list(next_current)

    def optimal_assignment(self, assignments, required_nodes=None):
        index_pair_sets = self.build_unique_index_pairs(assignments)
        best_score = -float('inf')
        best_mapping = None
        for assignment in index_pair_sets:
            current_score = 0
            has_nodes = set()
            for ix in assignment:
                current_score += assignments[ix]
                has_nodes.add(ix[0])
            if current_score > best_score and has_nodes == required_nodes:
                best_score = current_score
                best_mapping = assignment
        return best_mapping, best_score


topological_inclusion = TopologicalInclusionMatcher.compare


def exact_ordering_inclusion(self, other, substituents=True, tolerance=0, visited=None):
    '''
    A generalization of :meth:`~glypy.structure.monosaccharide.Monosaccharide.exact_ordering_equality` which
    allows for ``self`` to be matched to ``other``, but for ``other`` to include more. Consequently, this method is
    not commutative.
    '''
    if visited is None:
        visited = set()
    if (self.id, other.id) in visited:
        return True
    similar = commutative_similarity(self, other, tolerance, include_substituents=substituents)
    if similar:
        if substituents:
            other_substituents = dict(other.substituents())
            for a_pos, a_sub in self.substituents():
                b_sub = other_substituents.get(a_pos)
                if b_sub is None:  # pragma: no cover
                    return False
                if a_sub != b_sub:  # pragma: no cover
                    return False
        other_mods = dict(other.modifications.items())
        for a_pos, a_mod in self.modifications.items():
            b_mod = other_mods.get(a_pos)
            if b_mod is None:  # pragma: no cover
                return False
            if a_mod != b_mod:  # pragma: no cover
                return False
        other_children = dict(other.children())
        for pos, a_child in self.children():
            b_child = other_children.get(pos)
            if b_child is None:  # pragma: no cover
                return False
            if a_child[0] == b_child[0]:
                if not exact_ordering_inclusion(a_child, b_child, substituents=substituents,
                                                tolerance=tolerance, visited=visited):
                    return False
            else:
                return False
        return True


def subtree_of(subtree, tree, exact=False, tolerance=0):
    '''
    Test to see if `subtree` is included in `tree` anywhere. Returns the
    node id number of the first occurence of `subtree` included in `tree` or |None|
    if it is not found.

    Parameters
    ----------
    subtree: Glycan
        The structure to search for. The search attempts to match the complete structure of subtree.
    tree: Glycan
        The sturcture to search in. The search iterates over each residue in `tree` and calls a comparator
        function, comparing the `subtree` to the substructure rooted at that residue.
    exact: bool
        If |True|, use :func:`exact_ordering_inclusion` to compare nodes. Otherwise use :func:`topological_inclusion`.
        Defaults to |False|.

    Returns
    -------
    |int| or |None| if no match

    '''
    if exact:
        comparator = exact_ordering_inclusion
    else:
        comparator = topological_inclusion
    tree_root = root(subtree)
    for node in tree:
        if comparator(tree_root, node):
            return node.id
    return None


def find_matching_subtree_roots(subtree, tree, exact=False):
    if exact:
        comparator = exact_ordering_inclusion
    else:
        comparator = topological_inclusion
    tree_root = root(subtree)
    matched_nodes = []
    for node in tree:
        if comparator(tree_root, node):
            matched_nodes.append(node)
    return matched_nodes


def walk_with(query, reference, visited=None, comparator=commutative_similarity):
    if visited is None:
        visited = set()
    query_root = root(query)
    reference_root = root(reference)
    node_stack = deque([[reference_root, query_root]])

    while len(node_stack) != 0:
        rnode, qnode = node_stack.pop()
        key = (rnode.id, qnode.id)
        if key in visited:
            continue
        visited.add(key)
        yield rnode, qnode

        for p, rlink in rnode.links.items():
            rparent = rlink.parent
            qparent = None
            rchild = rlink.child
            qchild = None
            if p != -1:
                for qlink in qnode.links[p]:
                    if comparator(qlink.parent, rparent):
                        qparent = qlink.parent
                        break
                if qparent is not None:
                    key = (rparent.id, qparent.id)
                    if key not in visited:
                        node_stack.append((rparent, qparent))
                for qlink in qnode.links[p]:
                    if comparator(qlink.child, rchild):
                        qchild = qlink.child
                        break
                if qchild is not None:
                    key = (rchild.id, qchild.id)
                    if key not in visited:
                        node_stack.append((rchild, qchild))
            else:
                for qlink in qnode.links.values():
                    if comparator(qlink.parent, rparent) and comparator(qlink.child, rchild):
                        if rnode is rparent:
                            if (rchild.id, qlink.child.id) in visited:
                                continue
                        elif rnode is rchild:
                            if (rparent.id, qlink.parent.id) in visited:
                                continue
                        qparent = qlink.parent
                        qchild = qlink.child
                        break
                if qparent is not None:
                    key = (rparent.id, qparent.id)
                    if key not in visited:
                        node_stack.append((rparent, qparent))
                    key = (rchild.id, qchild.id)
                    if key not in visited:
                        node_stack.append((rchild, qchild))
