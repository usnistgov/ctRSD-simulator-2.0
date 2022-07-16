

Functions
=========

This section lists the different functions available within ctRSD-simulator-2.0.


.. admonition:: Note:

   Before all available functions can be accessed the model must be downloaded, imported, and instantiated!

   More information on these steps can be found :ref:`here <ImportSim>`.


.. _global_rate_constants:

global_rate_constants
---------------------

global_rate_constants(*krz='False'*, *krsd='False'*, *krev='False'*, *krep='False'*, *krepr='False'*, *kth='False'*, *krzTh='False'*
, *krsdF='False'*, *krevF='False'*, *krevA='False'*, *krsdA='False'*, *krzA='False'*, *krevCG='False'*, 
*krsdCG='False'*, *krzCG='False'*, *ktxnO='False'*, *ktxnG='False'*, *ktxnTh='False'* , *ktxnF='False'*, *ktxnAG='False'*
, *ktxnCG='False'*, *ktxn='False'*, *kssdO='False'*, *kssdF='False'*, *kdsduG='False'*, *kdsdG='False'*, *kdsdGO='False'*
, *kdsduAG='False'*, *kdsdAG='False*', *kdsduCG='False'*, *kdsdCG='False'*, *kdrd='False'*, *kdeg='False'*, *kssd='False'*, *kdsd='False'*)

global_rate_constants is used to globally change rate constants for all species of the same type (single matrix), or certain groups of different types (multiple matrices), instead of changing a specific rate constant for just one individual species (index of a matrix).

**Parameters:**
	krz: *float*, *optional* 
		Base ctRSD gate ribozymal cleavage rate.

	krsd: *float*, *optional* 
		Base ctRSD gate forward reaction rate.

	krev: *float*, *optional* 
		Base ctRSD output reverse reaction rate.

	krep: *float*, *optional* 
		Reporter reaction rate.

	krepr: *float*, *optional* 
		Reverse reporter reaction rate.

	kth: *float*, *optional* 
		Thresholding reaction rate.

	krzTh: *float*, *optional* 
		Thresholding ribozymal cleavage rate.

	krsdF: *float*, *optional* 
		Fuel forward reaction rate.

	krevF: *float*, *optional* 
		Fuel reverse reaction rate.

	krz: *float*, *optional* 
		Ribozymal cleavage rate.

	krevA: *float*, *optional* 
		ctRSD AND gate reverse reaction rate.

	krsdA: *float*, *optional* 
		ctRSD AND gate forward reaction rate.

	krzA: *float*, *optional* 
		ctRSD AND gate ribozymal cleavage rate.

	krevCG: *float*, *optional* 
		ctRSD comparator gate reverse reaction rate.

	krsdCG: *float*, *optional* 
		ctRSD comparator gate forward reaction rate.

	krzCG: *float*, *optional* 
		ctRSD comparator gate ribozymal cleavage rate.

	ktxnO: *float*, *optional* 
		Transcription rate for output.

	ktxnG: *float*, *optional* 
		Transcription rate for gate.

	ktxnTh: *float*, *optional* 
		Transcription rate for threshold.

	ktxnF: *float*, *optional* 
		Transcription rate for fuel.

	ktxnAG: *float*, *optional* 
		Transcription rate for AND gate.

	ktxnCG: *float*, *optional* 
		Transcription rate for comparator gate.

	ktxn: *float*, *optional*
		Globally change all transcription rates (ktxnO,ktxnG,ktxnTh,ktxnF,ktxnAG,ktxnCG)

	kssdO: *float*, *optional* 
		Degredation rate for output (single stranded).

	kssdF: *float*, *optional* 
		Degredation rate for fuel (single stranded).

	kdsduG: *float*, *optional* 
		Degredation rate for uncleaved gate (double stranded).

	kdsdG: *float*, *optional* 
		Degredation rate for gate (double stranded).

	kdsdGO: *float*, *optional* 
		Degredation rate for gate-output complex (double stranded).

	kdsduAG: *float*, *optional* 
		Degredation rate for uncleaved AND gate (double stranded).

	kdsdAG: *float*, *optional* 
		Degredation rate for AND gate (double stranded).

	kdsduCG: *float*, *optional* 
		Degredation rate for uncleaved comparator gate (double stranded).

	kdsdCG: *float*, *optional* 
		Degredation rate for comparator gate (double stranded).

	kdrd: *float*, *optional* 
		Degredation rate for RNA:DNA hybrids.

	kssd: *float*, *optional* 
		Degredation rate for all single stranded species (kssdO,kssdF).

	kdsd: *float*, *optional* 
		Degredation rate for all double stranded species (kdsduG,kdsdG,kdsdGO,kdsduAG,kdsdAG,kdsduCG,kdsdCG).

	kdeg: *float*, *optional* 
		Degredation rate for all species (kssdO,kssdF,kdsduG,kdsdG,kdsdGO,kdsduAG,kdsdAG,kdsduCG,kdsdCG,kdrd).







.. _molecular_species:

molecular_species
-----------------

molecular_species(*name*, *DNA_con=0*, *ic='False'*, *krz='False'*, *krsd='False'*, *krev='False'*, *krep='False'*, *krepr='False'*, *kth='False'*, *krzTh='False'*, *krsdF='False'*, *krevF='False'*, *krevA='False'*, *krsdA='False'*, *krzA='False'*, *krevCG='False'*, *krsdCG='False'*, *krzCG='False'*, *ktxnO='False'*, *ktxnG='False'*, *ktxnTh='False'*, *ktxnF='False'*, *ktxnAG='False'*, *ktxnCG='False'*, *kssdO='False'*, *kssdF='False'*, *kdsduG='False'*, *kdsdG='False'*, *kdsdGO='False'*, *kdsduAG='False'*, *kdsdAG='False'*, *kdsduCG='False'*, *kdsdCG='False'*, *kdrd='False'*)

molecular_species is used to initialize all species involved in the system being simulated.


.. admonition:: Warning!

   All optional rate constant inputs can only change corresponding species when inputted with those specific species. A warning message will be issued otherwise.

**Parameters:**
	name: *string*
		Name of species being initialized
			* Input -> I{domain} / IN{domain} / INP{domain} / INPUT{domain} (all options work, not case sensitive)
			* Gate -> G{domainI,domainO} / GATE{domainI,domainO} (all options work, not case sensitive)
			* Reporter -> R{domain}, REP{domain}, REPORTER{domain} (all options work, not case sensitive)
			* Output -> O{domainI,domainO} / OUT{domainI,domainO} / OUTPUT{domainI,domainO} (all options work, not case sensitive)
			* Uncleaved Gate -> uG{domainI,domainO} (not case sensitive)
			* Gate-Input Complex -> GI{domain} (not case sensitive)
			* Gate-Output Complex -> GO{domainI,domainO} (not case sensitive)
			* Reporter-Output Complex -> RO{domainI,domainO} (not case sensitive)
			* Output Reporter -> S{domain} (not case sensitive)
			* Uncleaved Threshold -> uTH{domain} (not case sensitive)
			* Threshold -> TH{domain} (not case sensitive)
			* Fuel -> F{domain} (not case sensitive)
			* Fuel Gate -> GF{domain} (not case sensitive)
			* Uncleaved AND Gate -> uAG{domainI,domainO} (not case sensitive)
			* AND Gate -> G{domainI1.domainI2,domainO} / GATE{domainI1.domainI2,domainO} / AG{domainI1.domainI2,domainO} (all options work, not case sensitive)
			* AND Gate-Output Complex A -> AGOa{domainI,domainO} (not case sensitive)
			* AND Gate-Output Complex B -> AGOb{domainI,domainO} (not case sensitive)
			* AND Gate Fuel Complex B -> AGF{domain} (not case sensitive)
			* Uncleaved Comparator Gate -> uCG{domainI,domainO} (not case sensitive)
			* Comparator Gate -> CG{domainI,domainO} (not case sensitive)
			* Comparator Gate-Output Complex A -> CGOa{domainI,domainO} (not case sensitive)
			* Comparator Gate-Output Complex B -> CGOb{domainI,domainO} (not case sensitive)

	DNA_con: *float*, *if NONE,default=0*
		DNA template concentration for inputed species. This and ic are the two ways a user can initialize a component being involved in the system. (Only applies to Input,Gate,Reporter,Output,Threshold,Fuel,GF,AG,CG)

	ic: *float*, *optional*
		Initial Concentration for inputted species. This and DNA_con are the two ways a user can initialize a component being involved in the system.

	krz: *float*, *optional* 
		Base ctRSD gate ribozymal cleavage rate.

	krsd: *float*, *optional* 
		Base ctRSD gate forward reaction rate.

	krev: *float*, *optional* 
		Base ctRSD output reverse reaction rate.

	krep: *float*, *optional* 
		Reporter reaction rate.

	krepr: *float*, *optional* 
		Reverse reporter reaction rate.

	kth: *float*, *optional* 
		Thresholding reaction rate.

	krzTh: *float*, *optional* 
		Thresholding ribozymal cleavage rate.

	krsdF: *float*, *optional* 
		Fuel forward reaction rate.

	krevF: *float*, *optional* 
		Fuel reverse reaction rate.

	krz: *float*, *optional* 
		Ribozymal cleavage rate.

	krevA: *float*, *optional* 
		ctRSD AND gate reverse reaction rate.

	krsdA: *float*, *optional* 
		ctRSD AND gate forward reaction rate.

	krzA: *float*, *optional* 
		ctRSD AND gate ribozymal cleavage rate.

	krevCG: *float*, *optional* 
		ctRSD comparator gate reverse reaction rate.

	krsdCG: *float*, *optional* 
		ctRSD comparator gate forward reaction rate.

	krzCG: *float*, *optional* 
		ctRSD comparator gate ribozymal cleavage rate.

	ktxnO: *float*, *optional* 
		Transcription rate for output.

	ktxnG: *float*, *optional* 
		Transcription rate for gate.

	ktxnTh: *float*, *optional* 
		Transcription rate for threshold.

	ktxnF: *float*, *optional* 
		Transcription rate for fuel.

	ktxnAG: *float*, *optional* 
		Transcription rate for AND gate.

	ktxnCG: *float*, *optional* 
		Transcription rate for comparator gate.

	kssdO: *float*, *optional* 
		Degredation rate for output (single stranded).

	kssdF: *float*, *optional* 
		Degredation rate for fuel (single stranded).

	kdsduG: *float*, *optional* 
		Degredation rate for uncleaved gate (double stranded).

	kdsdG: *float*, *optional* 
		Degredation rate for gate (double stranded).

	kdsdGO: *float*, *optional* 
		Degredation rate for gate-output complex (double stranded).

	kdsduAG: *float*, *optional* 
		Degredation rate for uncleaved AND gate (double stranded).

	kdsdAG: *float*, *optional* 
		Degredation rate for AND gate (double stranded).

	kdsduCG: *float*, *optional* 
		Degredation rate for uncleaved comparator gate (double stranded).

	kdsdCG: *float*, *optional* 
		Degredation rate for comparator gate (double stranded).

	kdrd: *float*, *optional* 
		Degredation rate for RNA:DNA hybrids.




.. _simulate: 

simulate
-----------------

simulate()

.. _output_concentration: 

output_concentration
--------------------

output_concentration()


.. _transcription_calibration: 

transcription_calibration
-------------------------

transcription_calibration()

