import sys, os, fileinput
import numpy as np

class GapPlot:
    def __init__(self, gap=0.0, threshold=None, bins=300):
        self.gap = gap
        self.threshold = threshold
        self.bins = bins

        self.states = list()
        self.transitions = list()

    def add_charge_state(self, charge, dfe):
        self.states.append( [charge, dfe] )

    def get_transition(self, state1, state2):
        #print state1, state2
        label = '(%i/%i)' % (state1[0], state2[0])
        trans = (state2[1] - state1[1]) / abs(state1[0] - state2[0])
        return label, trans

    def plot_levels(self):
        # Start by finding the lost energy charge state
        sorted_states = sorted(self.states, key=lambda dfe: dfe[1])

        step = self.gap / self.bins
        data = list()

        level = 0.0
        nsteps = 0
        energy = sorted_states[0][1]

        if (self.threshold is None):
            self.threshold = step

        for i in range(0, len(sorted_states)-1):
            state1 = sorted_states[i]
            state2 = sorted_states[i+1]
            transition = self.get_transition(state1, state2)
            #print transition, state1, state2
            if state2[0] > state1[0]:
                print 'The transition %s is losing electrons as Fermi level increases.' % transition[0]
                print 'This is unphysical, please check your DFE\'s and run again.'
                raise RuntimeError('Unphysical defect transition levels')

            switch = False
            while not switch:
#                if (abs(energy - state2[1]) <= self.threshold):
                if (abs(level - transition[1]) < self.threshold):
                    print 'Transition %s found at %4.2f eV.' % (transition[0], transition[1])
                    print 'DFE = %4.2f eV' % energy
                    print transition[1], energy
                    switch = True
                elif (level > self.gap):
                    switch = True # Level cannot be higher than the gap, so truncate it here
                else:
                    #print '%10.5f%10.5f%10.5f' % (energy - state2[1], level, energy)
                    data.append( [level, energy] )

                    energy += state1[0] * step
                    level += step
                    nsteps += 1

        final_state = sorted_states[-1]
        while (level <= self.gap):
            energy += final_state[0] * step
            data.append( [level, energy] )
            level += step

#       for datum in data:
#            print '%10.5f%10.5f' % (datum[0], datum[1])        

if __name__ == '__main__':
    datafile = sys.argv[1]
    data = open(datafile, 'r')

    gap = float(data.readline())
    gp = GapPlot(gap=float(gap), bins=300)

    eof = False
    while (not eof):
#        state = raw_input('Enter charge and DFE, or hit enter if done:\n')
        line = data.readline()
        if not line:
            eof = True
        else:
            state = line.split()
            gp.add_charge_state(int(state[0]), float(state[1]))

    gp.plot_levels()
