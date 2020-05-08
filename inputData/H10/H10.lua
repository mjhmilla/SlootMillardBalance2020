return {
metadata = {
	 {scaling_used = {"deLeva1996_segmentedTrunk"},
	 subject_age = {25.0},
	 subject_height = {1.75},
	 subject_weight = {79.00},
	 subject_gender = {"female"},
	 subject_pelvisWidth = {0.2240},
	 subject_hipCenterWidth = {0.1700},
	 subject_shoulderCenterWidth = {0.2600},
	 subject_heelAnkleXOffset = {0.0630},
	 subject_heelAnkleZOffset = {0.0790},
	 subject_shoulderNeckZOffset = {0.0450},
	 subject_footWidth = {0.1150},
	},
},
gravity = { 0, 0, -9.81,},
configuration = {
	axis_front = { 1, 0, 0,},
	axis_right = { 0, -1, 0,},
	axis_up = { 0, 0, 1,},
},
points = {
	{name = "Thigh_Proximal_R", body = "Thigh_R", point = {-0.074337, 0.000000, 0.037169,},},
	{name = "Thigh_Distal_R", body = "Thigh_R", point = {-0.074337, 0.000000, -0.092921,},},
	{name = "Heel_Medial_R", body = "Foot_R", point = {-0.063000, 0.040250, -0.079000,},},
	{name = "Heel_Lateral_R", body = "Foot_R", point = {-0.063000, -0.040250, -0.079000,},},
	{name = "ForeFoot_Medial_R", body = "Foot_R", point = {0.201815, 0.051750, -0.079000,},},
	{name = "ForeFoot_Lateral_R", body = "Foot_R", point = {0.201815, -0.051750, -0.079000,},},
	{name = "Thigh_Proximal_L", body = "Thigh_L", point = {-0.074337, 0.000000, 0.037169,},},
	{name = "Thigh_Distal_L", body = "Thigh_L", point = {-0.074337, 0.000000, -0.092921,},},
	{name = "Heel_Medial_L", body = "Foot_L", point = {-0.063000, -0.040250, -0.079000,},},
	{name = "Heel_Lateral_L", body = "Foot_L", point = {-0.063000, 0.040250, -0.079000,},},
	{name = "ForeFoot_Medial_L", body = "Foot_L", point = {0.201815, -0.051750, -0.079000,},},
	{name = "ForeFoot_Lateral_L", body = "Foot_L", point = {0.201815, 0.051750, -0.079000,},},
	{name = "DistalMetacarpal_Medial_R", body = "Hand_R", point = {-0.034314, 0.008579, -0.085785,},},
	{name = "DistalMetacarpal_Lateral_R", body = "Hand_R", point = {0.034314, 0.008579, -0.085785,},},
	{name = "DistalMetacarpal_Medial_L", body = "Hand_L", point = {-0.034314, -0.008579, -0.085785,},},
	{name = "DistalMetacarpal_Lateral_L", body = "Hand_L", point = {0.034314, -0.008579, -0.085785,},},
},
constraint_sets = {
},
frames = {
	{name   = "Pelvis",
	parent = "ROOT",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 9.851300,
		com = 	{ 0.000000, 0.000000, 0.092999,},
		inertia = 
			{{ 0.061901, 0.000000, 0.000000,},
			{ 0.000000, 0.053355, 0.000000,},
			{ 0.000000, 0.000000, 0.065086,},},
	},
	joint = 
		{{ 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000,},
		{ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000,},
		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_IAS = 	{ 0.137302, 0.091535, 0.109841,},
		R_IAS = 	{ 0.137302, -0.091535, 0.109841,},
		L_IPS = 	{ -0.128148, 0.054921, 0.128148,},
		R_IPS = 	{ -0.128148, -0.054921, 0.128148,},},
	visuals = {{
		src         = "pelvis.obj",
		dimensions  = 	{ 0.201600, 0.224000, 0.274604,},
		mesh_center = 	{ 0.000000, 0.000000, 0.045767,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Thigh_R",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, -0.085000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.676200,
		com = 	{ 0.000000, 0.000000, -0.134253,},
		inertia = 
			{{ 0.219637, 0.000000, 0.000000,},
			{ 0.000000, 0.213726, 0.000000,},
			{ 0.000000, 0.000000, 0.042333,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FTC = 	{ 0.000000, -0.074337, -0.122656,},
		R_FLE = 	{ 0.000000, -0.055753, -0.371686,},
		R_FME = 	{ 0.000000, 0.055753, -0.371686,},},
	visuals = {{
		src         = "thighR.obj",
		dimensions  = 	{ 0.148674, 0.148674, 0.446023,},
		mesh_center = 	{ 0.000000, 0.000000, -0.185843,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Shank_R",
	parent = "Thigh_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.371686,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.799900,
		com = 	{ 0.000000, 0.000000, -0.192529,},
		inertia = 
			{{ 0.053016, 0.000000, 0.000000,},
			{ 0.000000, 0.051440, 0.000000,},
			{ 0.000000, 0.000000, 0.006295,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FAX = 	{ 0.000000, -0.088478, -0.022120,},
		R_TTC = 	{ 0.022120, 0.000000, -0.066359,},
		R_FAL = 	{ 0.000000, -0.066359, -0.442392,},
		R_TAM = 	{ 0.000000, 0.066359, -0.442392,},},
	visuals = {{
		src         = "shankR.obj",
		dimensions  = 	{ 0.132718, 0.132718, 0.530870,},
		mesh_center = 	{ 0.000000, 0.000000, -0.221196,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Foot_R",
	parent = "Shank_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.442392,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.019100,
		com = 	{ 0.029432, 0.000000, -0.039500,},
		inertia = 
			{{ 0.004831, 0.000000, 0.000000,},
			{ 0.000000, 0.004206, 0.000000,},
			{ 0.000000, 0.000000, 0.001044,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_FM1 = 	{ 0.161192, 0.069082, -0.069082,},
		R_FM2 = 	{ 0.230274, 0.000000, -0.057568,},
		R_FM5 = 	{ 0.161192, -0.069082, -0.069082,},
		R_FCC = 	{ -0.069082, 0.000000, -0.023027,},},
	visuals = {{
		src         = "footR.obj",
		dimensions  = 	{ 0.299356, 0.115000, 0.059250,},
		mesh_center = 	{ 0.086678, 0.000000, -0.049375,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Thigh_L",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.085000, 0.000000,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.676200,
		com = 	{ 0.000000, 0.000000, -0.134253,},
		inertia = 
			{{ 0.219637, 0.000000, 0.000000,},
			{ 0.000000, 0.213726, 0.000000,},
			{ 0.000000, 0.000000, 0.042333,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FTC = 	{ 0.000000, 0.074337, -0.122656,},
		L_FLE = 	{ 0.000000, 0.055753, -0.371686,},
		L_FME = 	{ 0.000000, -0.055753, -0.371686,},},
	visuals = {{
		src         = "thighL.obj",
		dimensions  = 	{ 0.148674, 0.148674, 0.446023,},
		mesh_center = 	{ 0.000000, 0.000000, -0.185843,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Shank_L",
	parent = "Thigh_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.371686,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 3.799900,
		com = 	{ 0.000000, 0.000000, -0.192529,},
		inertia = 
			{{ 0.053016, 0.000000, 0.000000,},
			{ 0.000000, 0.051440, 0.000000,},
			{ 0.000000, 0.000000, 0.006295,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FAX = 	{ 0.000000, 0.088478, -0.022120,},
		L_TTC = 	{ 0.022120, 0.000000, -0.066359,},
		L_FAL = 	{ 0.000000, 0.066359, -0.442392,},
		L_TAM = 	{ 0.000000, -0.066359, -0.442392,},},
	visuals = {{
		src         = "shankL.obj",
		dimensions  = 	{ 0.132718, 0.132718, 0.530870,},
		mesh_center = 	{ 0.000000, 0.000000, -0.221196,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Foot_L",
	parent = "Shank_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.442392,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.019100,
		com = 	{ 0.029432, 0.000000, -0.039500,},
		inertia = 
			{{ 0.004831, 0.000000, 0.000000,},
			{ 0.000000, 0.004206, 0.000000,},
			{ 0.000000, 0.000000, 0.001044,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_FM1 = 	{ 0.161192, -0.069082, -0.069082,},
		L_FM2 = 	{ 0.230274, 0.000000, -0.057568,},
		L_FM5 = 	{ 0.161192, 0.069082, -0.069082,},
		L_FCC = 	{ -0.069082, 0.000000, -0.023027,},},
	visuals = {{
		src         = "footL.obj",
		dimensions  = 	{ 0.299356, 0.115000, 0.059250,},
		mesh_center = 	{ 0.086678, 0.000000, -0.049375,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "MiddleTrunk",
	parent = "Pelvis",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.183069,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 11.573500,
		com = 	{ 0.000000, 0.000000, 0.113643,},
		inertia = 
			{{ 0.093046, 0.000000, 0.000000,},
			{ 0.000000, 0.062191, 0.000000,},
			{ 0.000000, 0.000000, 0.085470,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		SXS = 	{ 0.122500, 0.000000, 0.103537,},
		MAI = 	{ -0.122500, 0.000000, 0.103537,},
		LV1 = 	{ -0.122500, 0.000000, 0.000000,},
		LV3 = 	{ -0.122500, 0.000000, -0.051769,},},
	visuals = {{
		src         = "middleTrunk.obj",
		dimensions  = 	{ 0.245000, 0.245000, 0.310612,},
		mesh_center = 	{ 0.000000, 0.000000, 0.051769,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperTrunk",
	parent = "MiddleTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.207075,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 12.205500,
		com = 	{ 0.000000, 0.000000, 0.113836,},
		inertia = 
			{{ 0.140176, 0.000000, 0.000000,},
			{ 0.000000, 0.063645, 0.000000,},
			{ 0.000000, 0.000000, 0.130135,},},
	},
	joint = 
		{},
	markers = {
		CV7 = 	{ -0.098000, 0.000000, 0.183977,},
		SJN = 	{ 0.102900, 0.000000, 0.137983,},
		TV2 = 	{ -0.098000, 0.000000, 0.137983,},},
	visuals = {{
		src         = "upperTrunk.obj",
		dimensions  = 	{ 0.245000, 0.273000, 0.241470,},
		mesh_center = 	{ 0.000000, 0.000000, 0.091988,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "Head",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, 0.229971,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 5.277200,
		com = 	{ 0.000000, 0.000000, 0.126812,},
		inertia = 
			{{ 0.023417, 0.000000, 0.000000,},
			{ 0.000000, 0.027748, 0.000000,},
			{ 0.000000, 0.000000, 0.021721,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HEAD = 	{ 0.000000, 0.073742, 0.147484,},
		R_HEAD = 	{ 0.000000, -0.073742, 0.147484,},
		SGL    = 	{ 0.073742, 0.000000, 0.147484,},},
	visuals = {{
		src         = "head.obj",
		dimensions  = 	{ 0.196646, 0.196646, 0.258097,},
		mesh_center = 	{ 0.000000, 0.000000, 0.098323,},
		color       = 	{ 0.750000, 0.750000, 0.750000,},
		},},
	},
	{name   = "UpperArm_R",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, -0.130000, 0.184971,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.014500,
		com = 	{ 0.000000, 0.000000, -0.159661,},
		inertia = 
			{{ 0.011987, 0.000000, 0.000000,},
			{ 0.000000, 0.010485, 0.000000,},
			{ 0.000000, 0.000000, 0.003397,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_SAE = 	{ 0.000000, 0.000000, 0.027748,},
		R_HUM = 	{ 0.000000, -0.041622, -0.185911,},
		R_HLE = 	{ 0.000000, -0.041622, -0.277478,},},
	visuals = {{
		src         = "upperArmR.obj",
		dimensions  = 	{ 0.138739, 0.110991, 0.305226,},
		mesh_center = 	{ 0.000000, 0.000000, -0.138739,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "LowerArm_R",
	parent = "UpperArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.277478,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.090200,
		com = 	{ 0.000000, 0.000000, -0.121536,},
		inertia = 
			{{ 0.005278, 0.000000, 0.000000,},
			{ 0.000000, 0.005117, 0.000000,},
			{ 0.000000, 0.000000, 0.000685,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_USP = 	{ -0.053317, 0.000000, -0.266585,},
		R_RSP = 	{ 0.053317, 0.000000, -0.266585,},},
	visuals = {{
		src         = "lowerArmR.obj",
		dimensions  = 	{ 0.106634, 0.079976, 0.293244,},
		mesh_center = 	{ 0.000000, 0.000000, -0.133293,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "Hand_R",
	parent = "LowerArm_R",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.266585,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.442400,
		com = 	{ 0.000000, 0.000000, -0.058797,},
		inertia = 
			{{ 0.000775, 0.000000, 0.000000,},
			{ 0.000000, 0.000563, 0.000000,},
			{ 0.000000, 0.000000, 0.000309,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		R_HM2 = 	{ 0.000000, 0.000000, -0.137256,},},
	visuals = {{
		src         = "handR.obj",
		dimensions  = 	{ 0.120099, 0.051471, 0.171571,},
		mesh_center = 	{ 0.000000, 0.000000, -0.085785,},
		color       = 	{ 0.700000, 0.200000, 0.300000,},
		},},
	},
	{name   = "UpperArm_L",
	parent = "UpperTrunk",
	joint_frame = {
		r = 	{ 0.000000, 0.130000, 0.184971,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 2.014500,
		com = 	{ 0.000000, 0.000000, -0.159661,},
		inertia = 
			{{ 0.011987, 0.000000, 0.000000,},
			{ 0.000000, 0.010485, 0.000000,},
			{ 0.000000, 0.000000, 0.003397,},},
	},
	joint = 
		{{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},
		{ 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_SAE = 	{ 0.000000, 0.000000, 0.027748,},
		L_HUM = 	{ 0.000000, 0.041622, -0.138739,},
		L_HLE = 	{ 0.000000, 0.041622, -0.277478,},},
	visuals = {{
		src         = "upperArmL.obj",
		dimensions  = 	{ 0.138739, 0.110991, 0.305226,},
		mesh_center = 	{ 0.000000, -0.000000, -0.138739,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "LowerArm_L",
	parent = "UpperArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.277478,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 1.090200,
		com = 	{ 0.000000, 0.000000, -0.121536,},
		inertia = 
			{{ 0.005278, 0.000000, 0.000000,},
			{ 0.000000, 0.005117, 0.000000,},
			{ 0.000000, 0.000000, 0.000685,},},
	},
	joint = 
{		{ 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_USP = 	{ -0.053317, 0.000000, -0.266585,},
		L_RSP = 	{ 0.053317, 0.000000, -0.266585,},},
	visuals = {{
		src         = "lowerArmL.obj",
		dimensions  = 	{ 0.106634, 0.079976, 0.293244,},
		mesh_center = 	{ 0.000000, 0.000000, -0.133293,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	{name   = "Hand_L",
	parent = "LowerArm_L",
	joint_frame = {
		r = 	{ 0.000000, 0.000000, -0.266585,},
		E = 
			{{ 1.000000, 0.000000, 0.000000,},
			{ 0.000000, 1.000000, 0.000000,},
			{ 0.000000, 0.000000, 1.000000,},},
	},
	body = {
		mass   = 0.442400,
		com = 	{ 0.000000, 0.000000, -0.058797,},
		inertia = 
			{{ 0.000775, 0.000000, 0.000000,},
			{ 0.000000, 0.000563, 0.000000,},
			{ 0.000000, 0.000000, 0.000309,},},
	},
	joint = 
{		{ 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,},},
	markers = {
		L_HM2 = 	{ 0.000000, 0.000000, -0.137256,},},
	visuals = {{
		src         = "handL.obj",
		dimensions  = 	{ 0.120099, 0.051471, 0.171571,},
		mesh_center = 	{ 0.000000, 0.000000, -0.085785,},
		color       = 	{ 0.200000, 0.700000, 0.300000,},
		},},
	},
	},
}