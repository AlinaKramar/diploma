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


def generate(n_parents, n_children, cis_fraction=.5, n_allels=100):
    parents = generate_parents(n_parents, n_allels)
    children = cross(parents, n_children)
    orgamisms = parents + children
    random.shuffle(orgamisms)
    return orgamisms


def generate_parents(n, n_allels):
    def random_parent(id):
        allels = [(random.choice([1, 2]), random.choice([1, 2])) for _ in range(n_allels)]
        sex = random.choice([MALE, FEMALE])
        return Organism(id, sex, allels)

    return [random_parent(id) for id in range(1, n + 1)]


def cross(parents, n_children, start_id=None):
    if not start_id:
        start_id = max(x.id for x in parents) + 1
    mothers = [m for m in parents if m.sex == FEMALE]
    fathers = [f for f in parents if f.sex == MALE]
    assert mothers and fathers, "All parents have the same sex"

    def cross_pair(id):
        sex = random.choice([MALE, FEMALE])
        m = random.choice(mothers)
        f = random.choice(fathers)
        assert len(m.allels) == len(f.allels), "Incompatible parents"
        mi = random.choice([0, 1])
        fi = random.choice([0, 1])
        allels = []
        for (m_a, f_a) in zip(m.allels, f.allels):
            if random.random() < 0.2:
                mi = random.choice([0, 1])
                fi = random.choice([0, 1])
            allels.append((m_a[mi], f_a[fi]))

        return Organism(id, sex, allels, m, f)

    return [cross_pair(id) for id in range(start_id, start_id + n_children)]

def to_file(orgamisms, destination="out2.gen"):

    def get_header():
        n_allels = len(orgamisms[0].allels)
        return "1\n{n_allels}\n{allel_names}\n\n1\n{n_orgnisms}\n".format(
            n_allels=n_allels,
            allel_names="x "*n_allels,
            n_orgnisms=len(orgamisms)
        )

    def format_allels(allels):
        return " ".join("{} {}".format(a, b) for (a, b) in allels)


    with open(destination, "w") as f:
        f.write(get_header())
        for o in orgamisms:
            s = "{id} {mother} {father} {sex} \n{allels}\n".format(
                id=o.id,
                mother=o.mother.id if o.mother else 0,
                father=o.father.id if o.father else 0,
                sex={MALE: 0, FEMALE: 1}[o.sex],
                allels=format_allels(o.allels)
            )
            f.write(s)

if __name__=="__main__":
    to_file(generate(n_parents=3, n_children=7, n_allels=35))
