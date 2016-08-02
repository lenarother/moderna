#!/usr/bin/env python
"""
Allowed and disallowed geometry values in RNA.
"""

class GeometryStandards:
    """Defines allowed and disallowed geometry values."""
    bonds = {
            "X:P,X:OP1":[( 1.40762, 1.72316)],
            "X:P,X:OP2":[( 1.38361, 1.70171)],
            "X:P,X:O5'":[( 1.50116, 1.8041)], # 2.08041
            "X:O3',X+1:P":[( 1.50116, 1.8041)], # 2.08041
            "X:O5',X:C5'":[( 1.28865, 1.62451)],
            "X:C5',X:C4'":[( 1.38598, 1.60948)],
            "X:C4',X:O4'":[( 1.37066, 1.54065)],
            "X:C4',X:C3'":[( 1.41850, 1.66191)],
            "X:O4',X:C1'":[( 1.34674, 1.57400)],
            "X:C3',X:O3'":[( 1.24104,  1.62451)],
            "X:C3',X:C2'":[( 1.36441, 1.61662)],
            "X:C2',X:O2'":[( 1.27818, 1.68238)],
            "X:C2',X:C1'":[( 1.41341, 1.63283)],
            "X[A]:C1',X[A]:N9":[( 1.39697, 1.53180)],
            "X[A]:N9,X[A]:C8":[( 1.30771, 1.47345)],
            "X[A]:C8,X[A]:N7":[( 1.24335, 1.47325)],
            "X[A]:N7,X[A]:C5":[( 1.27894, 1.42486)],
            "X[A]:C6,X[A]:N6":[( 1.28686, 1.38031)],
            "X[A]:C6,X[A]:C5":[( 1.34945, 1.47031)],
            "X[A]:C5,X[A]:C4":[( 1.30790, 1.45810)],
            "X[A]:C4,X[A]:N9":[( 1.31125, 1.41176)],
            "X[A]:C6,X[A]:N1":[( 1.28960, 1.40506)],
            "X[A]:N1,X[A]:C2":[( 1.29281, 1.38839)],
            "X[A]:C2,X[A]:N3":[( 1.28138, 1.51326)],
            "X[A]:N3,X[A]:C4":[( 1.28027, 1.37878)],
            "X[G]:C1',X[G]:N9":[( 1.27969, 1.60490)],
            "X[G]:N9,X[G]:C8":[( 1.29044, 1.47864)],
            "X[G]:C8,X[G]:N7":[( 1.27993, 1.49210)],
            "X[G]:N7,X[G]:C5":[( 1.34130, 1.49071)],
            "X[G]:C6,X[G]:O6":[( 1.18068, 1.29281)],
            "X[G]:C6,X[G]:C5":[( 1.35673, 1.57448)],
            "X[G]:C5,X[G]:C4":[( 1.32139, 1.46771)],
            "X[G]:C4,X[G]:N9":[( 1.31662, 1.47820)],
            "X[G]:C6,X[G]:N1":[( 1.27060, 1.43139)],
            "X[G]:N1,X[G]:C2":[( 1.32242, 1.42653)],
            "X[G]:C2,X[G]:N3":[( 1.29541, 1.41449)],
            "X[G]:N3,X[G]:C4":[( 1.30063, 1.39016)],
            "X[G]:C2,X[G]:N2":[( 1.22763, 1.43022)],
            "X[C]:C6,X[C]:C5":[( 1.26580, 1.45753)],
            "X[C]:C2,X[C]:O2":[( 1.16101, 1.30727)],
            "X[C]:C4,X[C]:N4":[( 1.27346, 1.40568)],
            "X[C]:C5,X[C]:C4":[( 1.31144, 1.52704)],
            "X[C]:C6,X[C]:N1":[( 1.22361, 1.43148)],
            "X[C]:N3,X[C]:C4":[( 1.25618, 1.41959)],
            "X[C]:N1,X[C]:C2":[( 1.29249, 1.47254)],
            "X[C]:C2,X[C]:N3":[( 1.30341, 1.43501)],
            "X[U]:C6,X[U]:C5":[( 1.30756, 1.40029)],
            "X[U]:C2,X[U]:O2":[( 1.17294, 1.30953)],
            "X[U]:C4,X[U]:O4":[( 1.17314, 1.30663)],
            "X[U]:C5,X[U]:C4":[( 1.38859, 1.53447)],
            "X[U]:C6,X[U]:N1":[( 1.32523, 1.44609)],
            "X[U]:N3,X[U]:C4":[( 1.32431, 1.45215)],
            "X[U]:N1,X[U]:C2":[( 1.33628, 1.46292)],
            "X[U]:C2,X[U]:N3":[( 1.33838, 1.40959)],
    }
    angles = {
            "X:P,X:O5',X:C5'":[(99.59706,139.93934)],
            "X:O5',X:C5',X:C4'":[(89.73461,131.30779)],
            "X:C5',X:C4',X:C3'":[(98.52296,131.88271)],
            "X:C4',X:C3',X:O3'":[(88.86713,127.18103)],
            "X:C3',X:O3',X+1:P":[( 100.0,140.0)],
            "X:O3',X+1:P,X+1:O5'":[(90.0,120.0)],
            "X:C5',X:C4',X:O4'":[(97.85665,146.68981)],
            "X:C4',X:O4',X:C1'":[(84.21729,115.71136)],
            "X:O4',X:C4',X:C3'":[(94.63959,110.53526)],
            "X:C4',X:C3',X:C2'":[(94.17150,112.81071)],
            "X:C3',X:C2',X:O2'":[(80.03621,126.84823)],
            "X:O3',X:C3',X:C2'":[(77.00673,138.55925)],
            "X:C3',X:C2',X:C1'":[(80.81583,108.24105)],
            "X:O2',X:C2',X:C1'":[(97.87503,136.10792)],
            "X:C2',X:C1',X:O4'":[(99.24436,131.37196)],
            "X[A]:C1',X[A]:N9,X[A]:C8":[(121.29182,134.87297)],
            "X[A]:N9,X[A]:C8,X[A]:N7":[(100.11655,118.57393)],
            "X[A]:C8,X[A]:N7,X[A]:C5":[(101.58735,107.94879)],
            "X[A]:N7,X[A]:C5,X[A]:C6":[(126.96407,137.17149)],
            "X[A]:N7,X[A]:C5,X[A]:C4":[(107.82678,114.80066)],
            "X[A]:C5,X[A]:C4,X[A]:N9":[(102.28972,108.88976)],
            "X[A]:C5,X[A]:C4,X[A]:N3":[(121.52363,132.73857)],
            "X[A]:C5,X[A]:C6,X[A]:N1":[(112.41227,122.64358)],
            "X[A]:C6,X[A]:N1,X[A]:C2":[(115.13842,125.84871)],
            "X[A]:N1,X[A]:C2,X[A]:N3":[(122.06095,132.72824)],
            "X[A]:C2,X[A]:N3,X[A]:C4":[(103.78598,116.86889)],
            "X[A]:N3,X[A]:C4,X[A]:N9":[(122.36966,130.87604)],
            "X[A]:N6,X[A]:C6,X[A]:C5":[(116.58830,128.68632)],
            "X[A]:N6,X[A]:C6,X[A]:N1":[(114.43765,126.64537)],
            "X[G]:C1',X[G]:N9,X[G]:C8":[(112.01297,135.87611)],
            "X[G]:N9,X[G]:C8,X[G]:N7":[(97.88571,117.97305)],
            "X[G]:C8,X[G]:N7,X[G]:C5":[(97.52309,109.50209)],
            "X[G]:N7,X[G]:C5,X[G]:C6":[(126.78408,137.42652)],
            "X[G]:N7,X[G]:C5,X[G]:C4":[(108.17958,116.01254)],
            "X[G]:C5,X[G]:C4,X[G]:N9":[(98.24771,107.90632)],
            "X[G]:C5,X[G]:C4,X[G]:N3":[(123.79833,132.56276)],
            "X[G]:C5,X[G]:C6,X[G]:N1":[(108.46936,117.11570)],
            "X[G]:C6,X[G]:N1,X[G]:C2":[(119.75726,129.35773)],
            "X[G]:N1,X[G]:C2,X[G]:N3":[(119.36242,126.99662)],
            "X[G]:C2,X[G]:N3,X[G]:C4":[(108.84509,116.49279)],
            "X[G]:N3,X[G]:C4,X[G]:N9":[(122.86395,131.90886)],
            "X[G]:O6,X[G]:C6,X[G]:C5":[(117.66282,134.75018)],
            "X[G]:O6,X[G]:C6,X[G]:N1":[(111.26501,128.96772)],
            "X[G]:N2,X[G]:C2,X[G]:N1":[(108.22873,123.31994)],
            "X[G]:N2,X[G]:C2,X[G]:N3":[(115.24447,125.40315)],
            "X[C]:N1,X[C]:C2,X[C]:N3":[(115.02983,122.99670)],
            "X[C]:C2,X[C]:N3,X[C]:C4":[(116.62450,124.32591)],
            "X[C]:N3,X[C]:C4,X[C]:N4":[(109.43859,124.31410)],
            "X[C]:C5,X[C]:C4,X[C]:N3":[(110.19781,127.90766)],
            "X[C]:C5,X[C]:C6,X[C]:N1":[(116.83415,127.04146)],
            "X[C]:C6,X[C]:N1,X[C]:C2":[(110.72263,125.76954)],
            "X[C]:C4,X[C]:C5,X[C]:C6":[(108.37411,121.48183)],
            "X[C]:N4,X[C]:C4,X[C]:C5":[(114.39426,127.70705)],
            "X[C]:O2,X[C]:C2,X[C]:N1":[(112.80935,126.39081)],
            "X[C]:O2,X[C]:C2,X[C]:N3":[(115.54160,128.95443)],
            "X[U]:N1,X[U]:C2,X[U]:N3":[(109.42133,121.19854)],
            "X[U]:C2,X[U]:N3,X[U]:C4":[(118.83532,133.97496)],
            "X[U]:N3,X[U]:C4,X[U]:O4":[(110.82094,125.12576)],
            "X[U]:C5,X[U]:C4,X[U]:N3":[(111.98415,121.92467)],
            "X[U]:C5,X[U]:C6,X[U]:N1":[(118.66714,127.00344)],
            "X[U]:C6,X[U]:N1,X[U]:C2":[(117.54900,125.67797)],
            "X[U]:C4,X[U]:C5,X[U]:C6":[(115.23890,123.50324)],
            "X[U]:O4,X[U]:C4,X[U]:C5":[(118.57411,135.90212)],
            "X[U]:O2,X[U]:C2,X[U]:N1":[(111.88691,129.45480)],
            "X[U]:O2,X[U]:C2,X[U]:N3":[(113.84482,138.65998)],
    }
    dihedrals = {
            # values from Richardson 2008
            "X:P,X:O5',X:C5',X:C4'":[(50.0, 290.0)], # beta
            "X:O5',X:C5',X:C4',X:C3'":[( 20.0, 95.0), (140.0, 215.0), (260.0, 335.0)], # gamma
            "X:C5',X:C4',X:C3',X:O3'":[( 55.0, 110.0), (120.0, 175.0)], # delta
            "X:C4',X:C3',X:O3',X+1:P":[( 155.0, 310.0)], # epsilon
            "X:C3',X:O3',X+1:P,X+1:O5'":[( 25.0, 335.0)], # zeta
            "X:O3',X+1:P,X+1:O5',X+1:C5'":[( 25.0, 335.0)], # alpha
            # O2' angles to make sure chirality is right
            "X:O3',X:C3',X:C2',X:O2'":[(0.0, 70.0), (290.0, 360)],
            "X:O4',X:C1',X:C2',X:O2'":[(60.0, 200.0)],
            # chi angle: everything is allowed
            #"X[A]:C8,X[A]:N9,X[A]:C1',X[A]:C2'":[( 2.77662,359.86181)],
            #"X[G]:C8,X[G]:N9,X[G]:C1',X[G]:C2'":[( 3.01550,358.96777)],
            #"X[C]:C2,X[C]:N1,X[C]:C1',X[C]:C2'":[(31.27292,352.95195)],
            #"X[U]:C2,X[U]:N1,X[U]:C1',X[U]:C2'":[(34.95466,356.89562)],
            # should be redundant by flat angle constraints
            #"X:C5',X:C4',X:C3',X:C2'":[(158.80652,294.73936)],
            #"X:O3',X:C3',X:C2',X:C1'":[(63.47648,267.40341)],
            #"X:C5',X:C4',X:O4',X:C1'":[(81.03942,196.47102)],
    }
    
    def __init__(self):
        self.values = {}
        for b in self.bonds: self.values[b] = self.bonds[b]
        for a in self.angles: self.values[a] = self.angles[a]
        for d in self.dihedrals: self.values[d] = self.dihedrals[d]
        
    def is_outlier(self, descriptor, value):
        if descriptor not in self.values:
            return True
        for lower, upper in self.values[descriptor]:
            if lower <= value <= upper:
                return False
        return True
    
    def get_standard(self, descriptor):
        """Returns value between lower and upper limit."""
        lower, upper = self.values[descriptor][0]
        return (lower+upper) / 2.0



BACKBONE_DIST_MATRIX = { 
        ("P", "O5'"): GeometryStandards.bonds["X:P,X:O5'"][0], 
        ("O5'", "C5'"): GeometryStandards.bonds["X:O5',X:C5'"][0], 
        ("C5'", "C4'"): GeometryStandards.bonds["X:C5',X:C4'"][0], 
        ("C4'", "C3'"): GeometryStandards.bonds["X:C4',X:C3'"][0], 
        ("C3'", "O3'"): GeometryStandards.bonds["X:C3',X:O3'"][0], 
    }

PHOSPHATE_DIST_MATRIX = { 
        ("P", "OP1"): GeometryStandards.bonds["X:P,X:OP1"][0], 
        ("P", "OP2"): GeometryStandards.bonds["X:P,X:OP2"][0], 
    }

O3_P_DIST_LOW, O3_P_DIST_HI  = GeometryStandards.bonds["X:O3',X+1:P"][0]
