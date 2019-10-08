# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from itertools import cycle, product, combinations
from networkx import compose, has_path


class Balance:
    __slots__ = ()

    def clone_subgraphs(self, g):
        if not isinstance(g, CGRContainer):
            raise InvalidData('only CGRContainer acceptable')

        r_group = []
        x_group = {}
        r_group_clones = []
        newcomponents = []

        ''' search bond breaks and creations
        '''
        components, lost_bonds, term_atoms = self.__split_graph(g)
        lost_map = {x: y for x, y in lost_bonds}
        ''' extract subgraphs and sort by group type (R or X)
        '''
        x_terminals = set(lost_map.values())
        r_terminals = set(lost_map)

        for i in components:
            x_terminal_atom = x_terminals.intersection(i)
            if x_terminal_atom:
                x_group[x_terminal_atom.pop()] = i
                continue

            r_terminal_atom = r_terminals.intersection(i)
            if r_terminal_atom:
                r_group.append([r_terminal_atom, i])
                continue

            newcomponents.append(i)
        ''' search similar R groups and patch.
        '''
        tmp = g
        for i in newcomponents:
            for k, j in r_group:
                gm = GraphMatcher(j, i, node_match=self.__node_match_products,
                                  edge_match=self.__edge_match_products)
                ''' search for similar R-groups started from bond breaks.
                '''
                mapping = next((x for x in gm.subgraph_isomorphisms_iter() if k.issubset(x) and
                                all(x[y] in term_atoms for y in k)), None)
                if mapping:
                    r_group_clones.append([k, mapping])
                    tmp = compose(tmp, self.__remap_group(j, tmp, mapping)[0])
                    break

        ''' add lose X groups to R groups
        '''
        for i, j in r_group_clones:
            for k in i:
                remappedgroup, mapping = self.__remap_group(x_group[lost_map[k]], tmp, {})
                tmp = CGRcore.union(tmp, remappedgroup)
                tmp.add_edge(j[k], mapping[lost_map[k]], s_bond=1, sp_bond=(1, None))

        if r_group_clones:
            tmp.meta.update(g.meta)
            return tmp

        return tmp.copy()

    @classmethod
    def __split_graph(cls, g):
        g = g.copy()
        lost_bonds = []
        term_atoms = []

        for l, n, m in list(cls.__get_substitution_paths(g)):
            nl = (n, l)
            if nl in lost_bonds:
                continue

            lost_bonds.append(nl)
            g.remove_edge(n, l)
            g.remove_edge(n, m)

        for n, m in list(cls.__get_broken_paths(g)):
            if not any(has_path(g, *x) for x in product((y for x in lost_bonds for y in x), (n, m))):
                g.remove_edge(n, m)
                term_atoms.append(n)
                term_atoms.append(m)

        return CGRcore.split(g), lost_bonds, term_atoms

    @staticmethod
    def __get_substitution_paths(g):
        """
        get atoms paths from detached atom to attached

        :param g: CGRContainer
        :return: tuple of atoms numbers
        """
        for n, nbrdict in g.adjacency():
            for m, l in combinations(nbrdict, 2):
                nms = nbrdict[m]['sp_bond']
                nls = nbrdict[l]['sp_bond']
                if nms == (1, None) and nls == (None, 1):
                    yield m, n, l
                elif nms == (None, 1) and nls == (1, None):
                    yield l, n, m

    @staticmethod
    def __get_broken_paths(g):
        for m, n, attr in g.edges(data=True):
            if attr['sp_bond'] == (None, 1):
                yield n, m
