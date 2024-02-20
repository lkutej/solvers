{
	labelList cellBoundaryFaceCount(omega_.size(), 0);

//	dimensionedScalar nu1(transportProperties.lookup("nu"));
        Info<<"Viscosity "<<nu1<<endl;

	const fvPatchList& patches = mesh_.boundary();

	//- Initialise the near-wall epsilon field to zero

	forAll(patches, patchi)
	{
		const fvPatch& curPatch = patches[patchi];

		if (isType<wallFvPatch>(curPatch))

		{
			forAll(curPatch, facei)
			{
			label faceCelli = curPatch.faceCells()[facei];
	
			omega_[faceCelli] = 0.0;
			}
		}
	}

	forAll(patches, patchi)
	{
		const fvPatch& curPatch = patches[patchi];
		if (isType<wallFvPatch>(curPatch))
			{
				forAll(curPatch, facei)
				{
				label faceCelli = curPatch.faceCells()[facei];

				// For corner cells (with two boundary or more faces),
				// epsilon in the near-wall cell are calculated as an average

				cellBoundaryFaceCount[faceCelli]++;

				omega_[faceCelli] =

				2.0*nu1.value()/sqr(yr_[faceCelli]);

				// Info<< epsilon[faceCelli]<< endl;
				}
			}
		}

		// perform the averaging

		forAll(patches, patchi)
		{
			const fvPatch& curPatch = patches[patchi];
			if (isType<wallFvPatch>(curPatch))
			{
				forAll(curPatch, facei)
					{
						label faceCelli = curPatch.faceCells()[facei];
						omega_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
					}
			}
		}


		forAll(patches, patchi)
			{
				const fvPatch& p = patches[patchi];
				if (isType<wallFvPatch>(p))
				{
					omegaEqn().setValues
					(
						p.faceCells(),
						omega_.boundaryField()[patchi].patchInternalField()
					);
				}
			}

}

