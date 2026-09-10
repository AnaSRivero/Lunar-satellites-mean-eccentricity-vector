function propagator = getMoonPropagator(initialOrbit, degree, order, ThirdBody)
%GETMOONPROPAGATOR Build an Orekit numerical propagator around the Moon.
%
%   propagator = getMoonPropagator(initialOrbit, degree, order, ThirdBody)
%
%   initialOrbit - org.orekit.orbits.Orbit, expressed in a Moon-centered
%                  inertial frame (see example scripts).
%   degree, order - degree/order of the lunar gravity field (GRGM1200B,
%                  read from 'lunar_grail.gfc' in the orekit-data folder).
%   ThirdBody    - 1 to add Earth + Sun point-mass third-body attraction,
%                  0 to propagate under the lunar field alone.
%
%   Requires setup_orekit(...) to have been called beforehand in this
%   MATLAB session/worker.

    import org.orekit.bodies.CelestialBodyFactory;
    import org.orekit.propagation.numerical.NumericalPropagator;
    import org.orekit.propagation.SpacecraftState;
    import org.orekit.orbits.OrbitType;
    import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
    import org.orekit.forces.gravity.ThirdBodyAttraction;
    import org.orekit.forces.gravity.potential.GravityFieldFactory;
    import org.orekit.forces.gravity.potential.ICGEMFormatReader;

    % --- Central body: the Moon ---
    moon = CelestialBodyFactory.getMoon();
    MoonPA = moon.getBodyOrientedFrame();
    satMass = 400.0; % arbitrary: no drag/SRP is modeled, so mass does not affect the dynamics

    % --- Integrator (Dormand-Prince 853) ---
    minStep = 0.001;
    maxStep = 1000.0;
    positionTolerance = 1e-6;

    tolerances = NumericalPropagator.tolerances(positionTolerance, initialOrbit, OrbitType.CARTESIAN);

    import org.hipparchus.ode.nonstiff.DormandPrince853Integrator
    integrator = DormandPrince853Integrator(minStep, maxStep, tolerances(1), tolerances(2));

    % --- Propagator ---
    propagator = NumericalPropagator(integrator);
    propagator.setOrbitType(OrbitType.CARTESIAN);

    % --- Lunar gravity field (GRGM1200B, truncated to degree/order) ---
    GravityFieldFactory.clearPotentialCoefficientsReaders();
    moonReader = ICGEMFormatReader('lunar_grail.gfc', false);
    GravityFieldFactory.addPotentialCoefficientsReader(moonReader);
    moonProvider = GravityFieldFactory.getNormalizedProvider(degree, order);
    propagator.addForceModel(HolmesFeatherstoneAttractionModel(MoonPA, moonProvider));

    % --- Third-body perturbations (Earth + Sun) ---
    if ThirdBody == 1
        propagator.addForceModel(ThirdBodyAttraction(CelestialBodyFactory.getEarth()));
        propagator.addForceModel(ThirdBodyAttraction(CelestialBodyFactory.getSun()));
    end

    % --- Initial state ---
    propagator.setInitialState(SpacecraftState(initialOrbit, satMass));

end
