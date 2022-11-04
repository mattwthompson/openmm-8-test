import openmm
import openmm.app
import openmm.unit


minimize = True

for platform_index in range(openmm.Platform.getNumPlatforms()):

    pdb_file = openmm.app.PDBFile("topology.pdb")
    system = openmm.XmlSerializer.deserialize(open("system.xml", "r").read())
    integrator = openmm.VerletIntegrator(1.0 * openmm.unit.femtoseconds)
    platform = openmm.Platform.getPlatform(platform_index)

    print(f"Using platform '{platform.getName()}'")

    try:
        simulation = openmm.app.Simulation(
            topology=pdb_file.topology,
            system=system,
            integrator=integrator,
            platform=platform,
        )
    except openmm.OpenMMException as exception:
        if "No compatible" in str(exception):
            print(str(exception))
            continue
        else:
            raise exception

    simulation.context.setPositions(pdb_file.positions)
    simulation.context.setPeriodicBoxVectors(
        *pdb_file.topology.getPeriodicBoxVectors(),
    )

    for index, force in enumerate(system.getForces()):
        force.setForceGroup(index)

    for force_index in range(system.getNumForces()):
        state = simulation.context.getState(
            getEnergy=True,
            groups={force_index},
        )
        print(force_index, state.getPotentialEnergy())
        del state

    if minimize:
        simulation.minimizeEnergy()

        state = simulation.context.getState(
            getPositions=True,
        )

        with open(f"minimized_{platform.getName()}.pdb", "w") as f:
            openmm.app.PDBFile.writeFile(
                topology=pdb_file.topology,
                positions=state.getPositions(),
                file=f,
            )
