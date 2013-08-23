import sys, os, fileinput
import numpy as np

from MAST.utility import MASTError

class GapPlot:
    def __init__(self, gap=0.0, threshold=None, bins=300, dfe=None):
        self.gap = gap
        self.threshold = threshold
        self.bins = bins
        self.dfe = dfe

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
        step = self.gap / self.bins

        if (self.threshold is None):
            self.threshold = step

        for condition in self.dfe:
            print 'Plotting levels for %s conditions' % condition
            for defect in self.dfe[condition]:
                print 'Analyzing leves for defect %s.' % defect
                state = self.dfe[condition][defect]
                data = self.plot_defect(state, step)

                print 'Writing energies and grid to file %s-%s.txt' % (defect, condition.replace(' ', '_'))
                with open('%s-%s.txt' % (defect, condition.replace(' ', '_')), 'w') as datafile:
                    for datum in data:
                        datafile.write('%10.5f%10.5f\n' % (datum[0], datum[1]))

                print str()

    def plot_defect(self, state, step):
        # Start by finding the lost energy charge state
        sorted_states = sorted(state, key=lambda dfe: dfe[1])

        data = list()
        level = 0.0
        nsteps = 0
        energy = sorted_states[0][1]

        for i in range(0, len(sorted_states)-1):
            state1 = sorted_states[i]
            state2 = sorted_states[i+1]
            transition = self.get_transition(state1, state2)
            #print transition, state1, state2
            if state2[0] > state1[0]:
                print 'The transition %s at %.2f is losing electrons as Fermi level increases.' % (transition[0], transition[1])
                #print 'This is unphysical, please check your DFE\'s and run again.'
                print 'Skipping this transition'
#                raise MASTError(self.__class__.__name__, 'Unphysical defect transition levels')
                # Since state2 is higher in energy, we neglect here, and continue on from state1
                final_state = state1
            else:
                final_state = state2
                switch = False
                while not switch:
    #                if (abs(energy - state2[1]) <= self.threshold):
                    if (abs(level - transition[1]) < self.threshold):
                        print 'Transition %s found at %4.2f eV.' % (transition[0], transition[1])
                        #print 'DFE = %4.2f eV' % energy
                        #print transition[1], energy
                        switch = True
                    elif (level > self.gap):
                        switch = True # Level cannot be higher than the gap, so truncate it here
                    else:
                        #print '%10.5f%10.5f%10.5f' % (energy - state2[1], level, energy)
                        data.append( [level, energy] )

                        energy += state1[0] * step
                        level += step
                        nsteps += 1

        #final_state = sorted_states[-1]
        while (level <= self.gap):
            energy += final_state[0] * step
            data.append( [level, energy] )
            level += step

        return data
