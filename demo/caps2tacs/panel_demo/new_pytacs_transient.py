""""
new pytacs transient function if tacs is messed up
"""

def createBDFtransientProb(self, tInit=0.0, tFinal=35, numSteps=10, options={}, amplitude=None):
        """
        Automatically define tacs problem classes with loads using information contained in BDF file.
        This function assumes all loads are specified in the BDF and allows users to
        skip setting loads in Python.

        Returns
        ----------
        structProblems : dict[TACSProblem]
            Dictionary containing a predfined TACSProblem for every loadcase found int the BDF.
            The dictionary keys are the loadcase IDs from the BDF.

        Notes
        -----
        Currently only supports LOAD, FORCE, MOMENT, GRAV, PLOAD2, and PLOAD4 cards.
        Currently only supports staticProblem (SOL 101) and modalProblems (SOL 103)
        """

        if self.assembler is None:
            raise self._initializeError()

        # Make sure cross-referencing is turned on in pynastran
        if self.bdfInfo.is_xrefed is False:
            self.bdfInfo.cross_reference()
            self.bdfInfo.is_xrefed = True

        vpn = self.varsPerNode
        loads = self.bdfInfo.loads
        nloads = len(loads)

        structProblems = {}

        # If subcases have been added in Nastran, then subCase 0 should not be run
        if len(self.bdfInfo.subcases) > 1:
            skipCaseZero = True
        else:
            skipCaseZero = False

        # Loop through every load set and create a corresponding structural problem
        for subCase in self.bdfInfo.subcases.values():
            if skipCaseZero and subCase.id == 0:
                continue

            if 'SUBTITLE' in subCase.params:
                name = subCase.params['SUBTITLE'][0]
            else:
                name = 'load_set_%.3d' % (subCase.id)

            problem = self.createTransientProblem(name, tInit=tInit, tFinal=tFinal, numSteps=numSteps, options=options)
            timeSteps = problem.getTimeSteps()

            if 'LOAD' in subCase.params:
                loadsID = subCase.params['LOAD'][0]
                # Get loads and scalers for this load case ID
                loadSet, loadScale, _ = self.bdfInfo.get_reduced_loads(loadsID)
                # Loop through every load in set and add it to problem

                for loadInfo, orig_scale in zip(loadSet, loadScale):
                    # Add any point force or moment cards
                    
                    for step_i, time in enumerate(timeSteps):
                        if amplitude is not None:
                            scale = amplitude(time) * orig_scale
                        if loadInfo.type == 'FORCE' or loadInfo.type == 'MOMENT':
                            nodeID = loadInfo.node_ref.nid

                            loadArray = numpy.zeros(vpn)
                            if loadInfo.type == 'FORCE' and vpn >= 3:
                                loadArray[:3] += scale * loadInfo.scaled_vector
                            elif loadInfo.type == 'MOMENT' and vpn >= 6:
                                loadArray[3:6] += scale * loadInfo.scaled_vector
                            problem.addLoadToNodes(step_i, nodeID, loadArray, nastranOrdering=True)

                        # Add any gravity loads
                        elif loadInfo.type == 'GRAV':
                            inertiaVec = np.zeros(3, dtype=self.dtype)
                            inertiaVec[:3] = scale * loadInfo.scale * loadInfo.N
                            problem.addInertialLoad(step_i,inertiaVec)

                        # Add any pressure loads
                        # Pressure load card specific to shell elements
                        elif loadInfo.type == 'PLOAD2':
                            elemIDs = loadInfo.eids
                            pressure = scale * loadInfo.pressure
                            problem.addPressureToElements(step_i,elemIDs, pressure, nastranOrdering=True)

                        # Alternate more general pressure load type
                        elif loadInfo.type == 'PLOAD4':
                            self._addPressureFromPLOAD4(step_i,problem, loadInfo, scale)

                        else:
                            self._TACSWarning("Unsupported load type "
                                        f" '{loadInfo.type}' specified for load set number {loadInfo.sid}, skipping load")

            # append to list of structural problems
            structProblems[subCase.id] = problem

        return structProblems
