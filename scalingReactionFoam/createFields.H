    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar D
    (
        transportProperties.lookup("D")
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field psi, xi and chi\n" << endl;

    volScalarField xi // conserved component
    (
        IOobject
        (
            "xi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField chi // reactive component
    (
        IOobject
        (
            "chi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField xi2 // conserved component
    (
        IOobject
        (
            "xi2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField chi2 // reactive component
    (
        IOobject
        (
            "chi2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Reology coefficients
    const dictionary rheology(transportProperties.subDict("rheology"));

    dimensionedScalar nu
    (
        rheology.lookup("nu")
    );

    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    #include "createMRF.H"       //added
    #include "createFvOptions.H" //added

/////////////////////////////////////////////////////////////////
