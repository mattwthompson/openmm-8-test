import openmm
import openmm.app
import openmm.unit

for platform_index in range(openmm.Platform.getNumPlatforms()):

    pdb_file = openmm.app.PDBFile("topology.pdb")
    system = openmm.XmlSerializer.deserialize(open("system.xml", "r").read())
    integrator = openmm.VerletIntegrator(1.0 * openmm.unit.femtoseconds)
    platform = openmm.Platform.getPlatform(platform_index)

    print(f"Using platform '{platform.getName()}'")

    simulation = openmm.app.Simulation(
        topology=pdb_file.topology,
        system=system,
        integrator=integrator,
        platform=platform,
    )

    simulation.context.setPositions(pdb_file.positions)
    simulation.context.setPeriodicBoxVectors(*pdb_file.topology.getPeriodicBoxVectors())

    for force_index in range(system.getNumForces()):
        state = simulation.context.getState(
            getEnergy=True,
            groups={force_index},
        )
        print(force_index, state.getPotentialEnergy())
        del state
