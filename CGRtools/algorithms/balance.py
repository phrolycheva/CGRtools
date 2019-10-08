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
from itertools import chain
from functools import reduce
from operator import or_


class Balance:
    __slots__ = ()

    def balance_substituents(self):
        reactants = set(chain.from_iterable(self.reactants))
        products = set(chain.from_iterable(self.products))
        added = products - reactants
        removed = reactants - products
        common = reactants & products

        common_products = reduce(or_, self.products).substructure(common, as_query=True)
        # find added and removed groups
        removed_groups = []
        skinned = set()
        for m in self.reactants:
            a = removed.intersection(m)
            if a:
                mb = m._bonds
                sk = set()
                for x in a:
                    sk.update(common.intersection(mb[x]))
                if sk:
                    skinned.update(sk)
                    removed_groups.append(a)

        added_groups = []
        for m in self.products:
            a = added.intersection(m)
            if a:
                for s in m.substructure(a, as_query=True).split():
                    try:
                        mapping = next(s.get_mapping(common_products))
                    except StopIteration:
                        continue
                    else:
                        added_groups.append(mapping)

