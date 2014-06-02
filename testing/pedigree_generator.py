#!/usr/bin/env python3

import random
import textwrap

MALE = "male"
FEMALE = "female"

class Organism(object):
    def __init__(self, id, sex, allels=None, mother=None, father=None):
        self.id = id
        self.sex = sex
        self.mother = mother
        self.father = father
        self.allels = allels

    def __str__(self):
        return "Organism: id={id}".format(id=self.id)


def generate(n_parents, n_generations, n_children, cis_fraction=.5, n_allels=100):
    parents = generate_parents(n_parents, n_allels)
    children = cross(parents, n_generations, n_children)
    orgamisms = parents + children
    random.shuffle(orgamisms)
    return orgamisms


def generate_parents(n, n_allels):
    def random_parent(id):
        allels = [(random.choice([1, 2]), random.choice([1, 2]))
                  for _ in range(n_allels)]
        sex = random.choice([MALE, FEMALE])
        return Organism(id, sex, allels)

    return [random_parent(id) for id in range(1, n + 1)]


def cross(parents, n_generations, n_children):
    def cross_pair(id, m, f):
        assert len(m.allels) == len(f.allels), "Incompatible parents"
        sex = random.choice([MALE, FEMALE])
        mi = random.choice([0, 1])
        fi = random.choice([0, 1])
        allels = []
        for (m_a, f_a) in zip(m.allels, f.allels):
            if random.random() < 0.35:
                mi = random.choice([0, 1])
                fi = random.choice([0, 1])
            allels.append((m_a[mi], f_a[fi]))
            org = Organism(id, sex, allels, m, f)

        return org

    population = parents
    for gen in range(n_generations):
        start_id = max(x.id for x in population) + 1
        children_ids = range(start_id, start_id + n_children)
        mothers = [m for m in population if m.sex == FEMALE]
        fathers = [f for f in population if f.sex == MALE]
        assert mothers and fathers, "All parents have the same sex"
        children = [cross_pair(id,
                               random.choice(mothers),
                               random.choice(fathers))
                    for id in children_ids]
        population += children


    return population

def print_organisms(orgamisms):

    def get_header():
        n_allels = len(orgamisms[0].allels)
        return "1\n{n_allels}\n{allel_names}\n\n1\n{n_orgnisms}\n".format(
            n_allels=n_allels,
            allel_names=" ".join( "L{}".format(i) for i in range(n_allels) ),
            n_orgnisms=len(orgamisms)
        )

    def format_allels(allels):
        return " ".join("{} {}".format(a, b) for (a, b) in allels)

    print(get_header(), end='')
    for o in orgamisms:
        s = "{id} {mother} {father} {sex} \n{allels}\n".format(
            id=o.id,
            mother=o.mother.id if o.mother else 0,
            father=o.father.id if o.father else 0,
            sex={MALE: 0, FEMALE: 1}[o.sex],
            allels=format_allels(o.allels)
        )
        print(s, end='')


if __name__=="__main__":
    import sys
    n_organisms, n_markers = [int(x) for x in sys.argv[1:]]
    n_generations = 5
    n_children = int(n_organisms / n_generations)
    organisms = generate(
        n_parents=3, n_generations=n_generations,
        n_children=n_children, n_allels=n_markers)

    print_organisms(organisms)
