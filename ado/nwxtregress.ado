*! nwxtregress
capture program drop nwxtregress
program define nwxtregress, eclass
syntax varlist(ts min=2) [if] 	, 	///
		/// spatial weight matrix settings
		dvarlag(string)					/// name of spatial weights matrix for dependent variable
		[								/// now optional
		NOSParse 						/// do not use sparse matrix internally
		/// sampling options			
		draws(integer 2000)				/// number of griddy gibs draws
		gridlength(integer 10000)		/// grind length
		NOMIT(integer 500)				/// number of omitted draws
		BArrypace(numlist)				/// settings for BarryPace Trick, iterations, maxorder default: 50 100
		uselud							/// use LUD instead of Barry-Pace
		/// output
		direct 							/// show total, indirect and direct effects
		]

		if "`dvarlag'" != "" {
			_nwxtregsar `*'
		}


end
