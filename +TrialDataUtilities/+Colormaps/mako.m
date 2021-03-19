function cmap = mako(nColors)
% from Seaborn, (c) Michael Waskom.  (c.f. seaborn/seaborn/cm.py)

mako_lut = [ ...
    0.04503935, 0.01482344, 0.02092227; ...
    0.04933018, 0.01709292, 0.02535719; ...
    0.05356262, 0.01950702, 0.03018802; ...
    0.05774337, 0.02205989, 0.03545515; ...
    0.06188095, 0.02474764, 0.04115287; ...
    0.06598247, 0.0275665 , 0.04691409; ...
    0.07005374, 0.03051278, 0.05264306; ...
    0.07409947, 0.03358324, 0.05834631; ...
    0.07812339, 0.03677446, 0.06403249; ...
    0.08212852, 0.0400833 , 0.06970862; ...
    0.08611731, 0.04339148, 0.07538208; ...
    0.09009161, 0.04664706, 0.08105568; ...
    0.09405308, 0.04985685, 0.08673591; ...
    0.09800301, 0.05302279, 0.09242646; ...
    0.10194255, 0.05614641, 0.09813162; ...
    0.10587261, 0.05922941, 0.103854  ; ...
    0.1097942 , 0.06227277, 0.10959847; ...
    0.11370826, 0.06527747, 0.11536893; ...
    0.11761516, 0.06824548, 0.12116393; ...
    0.12151575, 0.07117741, 0.12698763; ...
    0.12541095, 0.07407363, 0.1328442 ; ...
    0.12930083, 0.07693611, 0.13873064; ...
    0.13317849, 0.07976988, 0.14465095; ...
    0.13701138, 0.08259683, 0.15060265; ...
    0.14079223, 0.08542126, 0.15659379; ...
    0.14452486, 0.08824175, 0.16262484; ...
    0.14820351, 0.09106304, 0.16869476; ...
    0.15183185, 0.09388372, 0.17480366; ...
    0.15540398, 0.09670855, 0.18094993; ...
    0.15892417, 0.09953561, 0.18713384; ...
    0.16238588, 0.10236998, 0.19335329; ...
    0.16579435, 0.10520905, 0.19960847; ...
    0.16914226, 0.10805832, 0.20589698; ...
    0.17243586, 0.11091443, 0.21221911; ...
    0.17566717, 0.11378321, 0.21857219; ...
    0.17884322, 0.11666074, 0.2249565 ; ...
    0.18195582, 0.11955283, 0.23136943; ...
    0.18501213, 0.12245547, 0.23781116; ...
    0.18800459, 0.12537395, 0.24427914; ...
    0.19093944, 0.1283047 , 0.25077369; ...
    0.19381092, 0.13125179, 0.25729255; ...
    0.19662307, 0.13421303, 0.26383543; ...
    0.19937337, 0.13719028, 0.27040111; ...
    0.20206187, 0.14018372, 0.27698891; ...
    0.20469116, 0.14319196, 0.28359861; ...
    0.20725547, 0.14621882, 0.29022775; ...
    0.20976258, 0.14925954, 0.29687795; ...
    0.21220409, 0.15231929, 0.30354703; ...
    0.21458611, 0.15539445, 0.31023563; ...
    0.21690827, 0.15848519, 0.31694355; ...
    0.21916481, 0.16159489, 0.32366939; ...
    0.2213631 , 0.16471913, 0.33041431; ...
    0.22349947, 0.1678599 , 0.33717781; ...
    0.2255714 , 0.1710185 , 0.34395925; ...
    0.22758415, 0.17419169, 0.35075983; ...
    0.22953569, 0.17738041, 0.35757941; ...
    0.23142077, 0.18058733, 0.3644173 ; ...
    0.2332454 , 0.18380872, 0.37127514; ...
    0.2350092 , 0.18704459, 0.3781528 ; ...
    0.23670785, 0.190297  , 0.38504973; ...
    0.23834119, 0.19356547, 0.39196711; ...
    0.23991189, 0.19684817, 0.39890581; ...
    0.24141903, 0.20014508, 0.4058667 ; ...
    0.24286214, 0.20345642, 0.4128484 ; ...
    0.24423453, 0.20678459, 0.41985299; ...
    0.24554109, 0.21012669, 0.42688124; ...
    0.2467815 , 0.21348266, 0.43393244; ...
    0.24795393, 0.21685249, 0.4410088 ; ...
    0.24905614, 0.22023618, 0.448113  ; ...
    0.25007383, 0.22365053, 0.45519562; ...
    0.25098926, 0.22710664, 0.46223892; ...
    0.25179696, 0.23060342, 0.46925447; ...
    0.25249346, 0.23414353, 0.47623196; ...
    0.25307401, 0.23772973, 0.48316271; ...
    0.25353152, 0.24136961, 0.49001976; ...
    0.25386167, 0.24506548, 0.49679407; ...
    0.25406082, 0.2488164 , 0.50348932; ...
    0.25412435, 0.25262843, 0.51007843; ...
    0.25404842, 0.25650743, 0.51653282; ...
    0.25383134, 0.26044852, 0.52286845; ...
    0.2534705 , 0.26446165, 0.52903422; ...
    0.25296722, 0.2685428 , 0.53503572; ...
    0.2523226 , 0.27269346, 0.54085315; ...
    0.25153974, 0.27691629, 0.54645752; ...
    0.25062402, 0.28120467, 0.55185939; ...
    0.24958205, 0.28556371, 0.55701246; ...
    0.24842386, 0.28998148, 0.56194601; ...
    0.24715928, 0.29446327, 0.56660884; ...
    0.24580099, 0.29899398, 0.57104399; ...
    0.24436202, 0.30357852, 0.57519929; ...
    0.24285591, 0.30819938, 0.57913247; ...
    0.24129828, 0.31286235, 0.58278615; ...
    0.23970131, 0.3175495 , 0.5862272 ; ...
    0.23807973, 0.32226344, 0.58941872; ...
    0.23644557, 0.32699241, 0.59240198; ...
    0.2348113 , 0.33173196, 0.59518282; ...
    0.23318874, 0.33648036, 0.59775543; ...
    0.2315855 , 0.34122763, 0.60016456; ...
    0.23001121, 0.34597357, 0.60240251; ...
    0.2284748 , 0.35071512, 0.6044784 ; ...
    0.22698081, 0.35544612, 0.60642528; ...
    0.22553305, 0.36016515, 0.60825252; ...
    0.22413977, 0.36487341, 0.60994938; ...
    0.22280246, 0.36956728, 0.61154118; ...
    0.22152555, 0.37424409, 0.61304472; ...
    0.22030752, 0.37890437, 0.61446646; ...
    0.2191538 , 0.38354668, 0.61581561; ...
    0.21806257, 0.38817169, 0.61709794; ...
    0.21703799, 0.39277882, 0.61831922; ...
    0.21607792, 0.39736958, 0.61948028; ...
    0.21518463, 0.40194196, 0.62059763; ...
    0.21435467, 0.40649717, 0.62167507; ...
    0.21358663, 0.41103579, 0.62271724; ...
    0.21288172, 0.41555771, 0.62373011; ...
    0.21223835, 0.42006355, 0.62471794; ...
    0.21165312, 0.42455441, 0.62568371; ...
    0.21112526, 0.42903064, 0.6266318 ; ...
    0.21065161, 0.43349321, 0.62756504; ...
    0.21023306, 0.43794288, 0.62848279; ...
    0.20985996, 0.44238227, 0.62938329; ...
    0.20951045, 0.44680966, 0.63030696; ...
    0.20916709, 0.45122981, 0.63124483; ...
    0.20882976, 0.45564335, 0.63219599; ...
    0.20849798, 0.46005094, 0.63315928; ...
    0.20817199, 0.46445309, 0.63413391; ...
    0.20785149, 0.46885041, 0.63511876; ...
    0.20753716, 0.47324327, 0.63611321; ...
    0.20722876, 0.47763224, 0.63711608; ...
    0.20692679, 0.48201774, 0.63812656; ...
    0.20663156, 0.48640018, 0.63914367; ...
    0.20634336, 0.49078002, 0.64016638; ...
    0.20606303, 0.49515755, 0.6411939 ; ...
    0.20578999, 0.49953341, 0.64222457; ...
    0.20552612, 0.50390766, 0.64325811; ...
    0.20527189, 0.50828072, 0.64429331; ...
    0.20502868, 0.51265277, 0.64532947; ...
    0.20479718, 0.51702417, 0.64636539; ...
    0.20457804, 0.52139527, 0.64739979; ...
    0.20437304, 0.52576622, 0.64843198; ...
    0.20418396, 0.53013715, 0.64946117; ...
    0.20401238, 0.53450825, 0.65048638; ...
    0.20385896, 0.53887991, 0.65150606; ...
    0.20372653, 0.54325208, 0.65251978; ...
    0.20361709, 0.5476249 , 0.6535266 ; ...
    0.20353258, 0.55199854, 0.65452542; ...
    0.20347472, 0.55637318, 0.655515  ; ...
    0.20344718, 0.56074869, 0.65649508; ...
    0.20345161, 0.56512531, 0.65746419; ...
    0.20349089, 0.56950304, 0.65842151; ...
    0.20356842, 0.57388184, 0.65936642; ...
    0.20368663, 0.57826181, 0.66029768; ...
    0.20384884, 0.58264293, 0.6612145 ; ...
    0.20405904, 0.58702506, 0.66211645; ...
    0.20431921, 0.59140842, 0.66300179; ...
    0.20463464, 0.59579264, 0.66387079; ...
    0.20500731, 0.60017798, 0.66472159; ...
    0.20544449, 0.60456387, 0.66555409; ...
    0.20596097, 0.60894927, 0.66636568; ...
    0.20654832, 0.61333521, 0.66715744; ...
    0.20721003, 0.61772167, 0.66792838; ...
    0.20795035, 0.62210845, 0.66867802; ...
    0.20877302, 0.62649546, 0.66940555; ...
    0.20968223, 0.63088252, 0.6701105 ; ...
    0.21068163, 0.63526951, 0.67079211; ...
    0.21177544, 0.63965621, 0.67145005; ...
    0.21298582, 0.64404072, 0.67208182; ...
    0.21430361, 0.64842404, 0.67268861; ...
    0.21572716, 0.65280655, 0.67326978; ...
    0.21726052, 0.65718791, 0.6738255 ; ...
    0.21890636, 0.66156803, 0.67435491; ...
    0.220668  , 0.66594665, 0.67485792; ...
    0.22255447, 0.67032297, 0.67533374; ...
    0.22458372, 0.67469531, 0.67578061; ...
    0.22673713, 0.67906542, 0.67620044; ...
    0.22901625, 0.6834332 , 0.67659251; ...
    0.23142316, 0.68779836, 0.67695703; ...
    0.23395924, 0.69216072, 0.67729378; ...
    0.23663857, 0.69651881, 0.67760151; ...
    0.23946645, 0.70087194, 0.67788018; ...
    0.24242624, 0.70522162, 0.67813088; ...
    0.24549008, 0.70957083, 0.67835215; ...
    0.24863372, 0.71392166, 0.67854868; ...
    0.25187832, 0.71827158, 0.67872193; ...
    0.25524083, 0.72261873, 0.67887024; ...
    0.25870947, 0.72696469, 0.67898912; ...
    0.26229238, 0.73130855, 0.67907645; ...
    0.26604085, 0.73564353, 0.67914062; ...
    0.26993099, 0.73997282, 0.67917264; ...
    0.27397488, 0.74429484, 0.67917096; ...
    0.27822463, 0.74860229, 0.67914468; ...
    0.28264201, 0.75290034, 0.67907959; ...
    0.2873016 , 0.75717817, 0.67899164; ...
    0.29215894, 0.76144162, 0.67886578; ...
    0.29729823, 0.76567816, 0.67871894; ...
    0.30268199, 0.76989232, 0.67853896; ...
    0.30835665, 0.77407636, 0.67833512; ...
    0.31435139, 0.77822478, 0.67811118; ...
    0.3206671 , 0.78233575, 0.67786729; ...
    0.32733158, 0.78640315, 0.67761027; ...
    0.33437168, 0.79042043, 0.67734882; ...
    0.34182112, 0.79437948, 0.67709394; ...
    0.34968889, 0.79827511, 0.67685638; ...
    0.35799244, 0.80210037, 0.67664969; ...
    0.36675371, 0.80584651, 0.67649539; ...
    0.3759816 , 0.80950627, 0.67641393; ...
    0.38566792, 0.81307432, 0.67642947; ...
    0.39579804, 0.81654592, 0.67656899; ...
    0.40634556, 0.81991799, 0.67686215; ...
    0.41730243, 0.82318339, 0.67735255; ...
    0.4285828 , 0.82635051, 0.6780564 ; ...
    0.44012728, 0.82942353, 0.67900049; ...
    0.45189421, 0.83240398, 0.68021733; ...
    0.46378379, 0.83530763, 0.6817062 ; ...
    0.47573199, 0.83814472, 0.68347352; ...
    0.48769865, 0.84092197, 0.68552698; ...
    0.49962354, 0.84365379, 0.68783929; ...
    0.5114027 , 0.8463718 , 0.69029789; ...
    0.52301693, 0.84908401, 0.69288545; ...
    0.53447549, 0.85179048, 0.69561066; ...
    0.54578602, 0.8544913 , 0.69848331; ...
    0.55695565, 0.85718723, 0.70150427; ...
    0.56798832, 0.85987893, 0.70468261; ...
    0.57888639, 0.86256715, 0.70802931; ...
    0.5896541 , 0.8652532 , 0.71154204; ...
    0.60028928, 0.86793835, 0.71523675; ...
    0.61079441, 0.87062438, 0.71910895; ...
    0.62116633, 0.87331311, 0.72317003; ...
    0.63140509, 0.87600675, 0.72741689; ...
    0.64150735, 0.87870746, 0.73185717; ...
    0.65147219, 0.8814179 , 0.73648495; ...
    0.66129632, 0.8841403 , 0.74130658; ...
    0.67097934, 0.88687758, 0.74631123; ...
    0.68051833, 0.88963189, 0.75150483; ...
    0.68991419, 0.89240612, 0.75687187; ...
    0.69916533, 0.89520211, 0.76241714; ...
    0.70827373, 0.89802257, 0.76812286; ...
    0.71723995, 0.90086891, 0.77399039; ...
    0.72606665, 0.90374337, 0.7800041 ; ...
    0.73475675, 0.90664718, 0.78615802; ...
    0.74331358, 0.90958151, 0.79244474; ...
    0.75174143, 0.91254787, 0.79884925; ...
    0.76004473, 0.91554656, 0.80536823; ...
    0.76827704, 0.91856549, 0.81196513; ...
    0.77647029, 0.921603  , 0.81855729; ...
    0.78462009, 0.92466151, 0.82514119; ...
    0.79273542, 0.92773848, 0.83172131; ...
    0.8008109 , 0.93083672, 0.83829355; ...
    0.80885107, 0.93395528, 0.84485982; ...
    0.81685878, 0.9370938 , 0.85142101; ...
    0.82483206, 0.94025378, 0.8579751 ; ...
    0.83277661, 0.94343371, 0.86452477; ...
    0.84069127, 0.94663473, 0.87106853; ...
    0.84857662, 0.9498573 , 0.8776059 ; ...
    0.8564431 , 0.95309792, 0.88414253; ...
    0.86429066, 0.95635719, 0.89067759; ...
    0.87218969, 0.95960708, 0.89725384 ];

if nargin == 0
    cmap = mako_lut;
else
    cmap = TrialDataUtilities.Color.evalColorMapAt(mako_lut, linspace(0, 1, nColors));
end