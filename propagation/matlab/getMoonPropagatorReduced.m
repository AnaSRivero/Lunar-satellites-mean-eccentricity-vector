function propagator = getMoonPropagatorReduced(initialOrbit, degree, order, extraTerms, orekitDataDir)
%GETMOONPROPAGATORREDUCED Numerical propagator with a reduced lunar gravity field.
%
%   propagator = getMoonPropagatorReduced(initialOrbit, degree, order, ...
%                    extraTerms, orekitDataDir)
%
%   Same integrator and settings as getMoonPropagator.m, but the lunar
%   gravity field contains only the harmonics retained in the analytical
%   solution (all zonal harmonics and the non-zonal harmonics of odd
%   degree, up to degree/order), optionally augmented with the even-degree
%   non-zonal harmonics listed in extraTerms. No third-body attraction is
%   included.
%
%   This is the reduced numerical model used in the paper to assess the
%   effect of the omitted even-degree non-zonal harmonics (Section 5).
%
%   initialOrbit  - org.orekit.orbits.Orbit in a Moon-centered inertial frame
%   degree, order - truncation of the lunar gravity field
%   extraTerms    - N-by-2 array of [degree, order] pairs whose C_lm and
%                   S_lm are added back, e.g. [2 2; 4 3]; zeros(0, 2) for none
%   orekitDataDir - folder with the orekit-data set and 'lunar_grail.gfc';
%                   the filtered gravity field file is written there
%                   (see buildSelectiveLunarGravityModel.m)
%
%   Requires setup_orekit(...) to have been called beforehand.

    import org.orekit.bodies.CelestialBodyFactory;
    import org.orekit.propagation.numerical.NumericalPropagator;
    import org.orekit.propagation.SpacecraftState;
    import org.orekit.orbits.OrbitType;
    import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
    import org.orekit.forces.gravity.potential.GravityFieldFactory;
    import org.orekit.forces.gravity.potential.ICGEMFormatReader;
    import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;

    % --- Central body: the Moon ---
    moon = CelestialBodyFactory.getMoon();
    MoonPA = moon.getBodyOrientedFrame();
    satMass = 400.0; % arbitrary: no drag/SRP is modeled, so mass does not affect the dynamics

    % --- Integrator (Dormand-Prince 853), as in getMoonPropagator.m ---
    minStep = 0.001;
    maxStep = 1000.0;
    positionTolerance = 1e-6;
    tolerances = NumericalPropagator.tolerances(positionTolerance, initialOrbit, OrbitType.CARTESIAN);
    integrator = DormandPrince853Integrator(minStep, maxStep, tolerances(1), tolerances(2));

    propagator = NumericalPropagator(integrator);
    propagator.setOrbitType(OrbitType.CARTESIAN);

    % --- Reduced lunar gravity field ---
    filteredModelName = buildSelectiveLunarGravityModel(orekitDataDir, degree, order, extraTerms);
    GravityFieldFactory.clearPotentialCoefficientsReaders();
    moonReader = ICGEMFormatReader(filteredModelName, false);
    GravityFieldFactory.addPotentialCoefficientsReader(moonReader);
    moonProvider = GravityFieldFactory.getNormalizedProvider(degree, order);
    propagator.addForceModel(HolmesFeatherstoneAttractionModel(MoonPA, moonProvider));

    % --- Initial state ---
    propagator.setInitialState(SpacecraftState(initialOrbit, satMass));

end
