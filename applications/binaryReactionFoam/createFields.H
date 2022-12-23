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

    // Diffusion coefficients
    dimensionedScalar Dpsi
    (
        transportProperties.lookup("Dpsi")
    );

    dimensionedScalar D
    (
        transportProperties.lookup("D")
    );

    // Reaction coefficients
    dimensionedScalar Kf
    (
        transportProperties.lookup("Kf")
    );


    // Reology coefficients
    const dictionary rheology(transportProperties.subDict("rheology"));

    dimensionedScalar nu0
    (
        rheology.lookup("nu")
    );

    dimensionedScalar cm
    (
        rheology.lookup("cm")
    );

    dimensionedScalar Ca
    (
        rheology.lookup("Ca")
    );

    dimensionedScalar Cb
    (
        rheology.lookup("Cb")
    );

    dimensionedScalar Cc
    (
        rheology.lookup("Cc")
    );

    // Stokes number
    dimensionedScalar St
    (
        transportProperties.lookup("St")
    );

    // scaling factor for source terms
    dimensionedScalar psifact
    (
        transportProperties.lookup("psifact")
    );


    dimensionedTensor I
    (
    "I",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    tensor(1,0,0,0,1,0,0,0,1)
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

    volScalarField psi // solid component
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

    Info<< "Creating viscosity field\n" << endl;

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu0 *
        (
            Ca
          + Cb*psi*pow( (scalar(1)- (psi/cm)), -scalar(1))
          + Cc*pow( (scalar(1)- (psi/cm)), -scalar(2))
        )
    );

    #include "createPhi.H"

    volVectorField Uslip
    (
        IOobject
        (
            "Uslip",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U*scalar(0)
    );

//    volVectorField Uslip("Uslip",U*0);
    surfaceScalarField phislip("phislip",phi*0);


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    #include "createMRF.H"       //added
    #include "createFvOptions.H" //added

/////////////////////////////////////////////////////////////////
