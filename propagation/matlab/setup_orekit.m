function setup_orekit(orekitJarDir, orekitDataDir)
%SETUP_OREKIT Add Orekit/Hipparchus jars to the MATLAB Java classpath and
%point Orekit's DataProvidersManager at a local orekit-data folder.
%
%   setup_orekit(orekitJarDir, orekitDataDir)
%
%   orekitJarDir  - folder containing orekit-*.jar and hipparchus-*.jar
%                   (developed with Orekit 13.1.x / Hipparchus 4.0.x; jar
%                   files are located by name pattern, so any matching
%                   version present in orekitJarDir is picked up)
%   orekitDataDir - folder containing the standard orekit-data physical
%                   data set (https://gitlab.orekit.org/orekit/orekit-data)
%                   plus the lunar gravity field file 'lunar_grail.gfc'
%                   (see ../data/README.md for where to get it)
%
%   Call this once per MATLAB session (and once per parfor worker, since
%   each worker has its own Java classpath / DataProvidersManager).

    cp = javaclasspath('-all');
    if ~any(contains(cp, 'orekit'))
        requiredJars = {'orekit-*.jar', 'hipparchus-core-*.jar', ...
            'hipparchus-geometry-*.jar', 'hipparchus-ode-*.jar', ...
            'hipparchus-fitting-*.jar', 'hipparchus-optim-*.jar', ...
            'hipparchus-stat-*.jar'};
        for k = 1:numel(requiredJars)
            found = dir(fullfile(orekitJarDir, requiredJars{k}));
            if isempty(found)
                error('setup_orekit:jarNotFound', ...
                    'No jar matching "%s" found in %s', requiredJars{k}, orekitJarDir);
            end
            javaaddpath(fullfile(found(1).folder, found(1).name));
        end
    end

    import org.orekit.data.DataContext
    import org.orekit.data.DirectoryCrawler

    DM = DataContext.getDefault().getDataProvidersManager();
    DM.clearProviders();
    DM.clearLoadedDataNames();
    if DM.getProviders().isEmpty()
        DM.addProvider(DirectoryCrawler(java.io.File(orekitDataDir)));
    end
end
