{
	labelList cellBoundaryFaceCount(epsilon_.size(), 0);

//	dimensionedScalar nu(transportProperties.lookup("nu"));

	const fvPatchList& patches = mesh_.boundary();

	dimensionedScalar nu1("nu1", dimensionSet(0,2,-1,0,0,0,0), 2.0e-05);

	//- Initialise the near-wall epsilon field to zero

	forAll(patches, patchi)
	{
		const fvPatch& curPatch = patches[patchi];

		if (isType<wallFvPatch>(curPatch))

		{
			forAll(curPatch, facei)
			{
			label faceCelli = curPatch.faceCells()[facei];
	
			f_[faceCelli] = 0.0;
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

				f_[faceCelli] =

				2.0*nu1.value()*zeta_[faceCelli]/sqr(yr_[faceCelli]);

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
						f_[faceCelli] /= cellBoundaryFaceCount[faceCelli];
					}
			}
		}


		forAll(patches, patchi)
			{
				const fvPatch& p = patches[patchi];
				if (isType<wallFvPatch>(p))
				{
					fEqn().setValues
					(
						p.faceCells(),
						f_.boundaryField()[patchi].patchInternalField()
					);
				}
			}

}

