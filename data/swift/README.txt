March 2nd, 2020
My understanding is that this delay depends purely on the conditions of the telemetry, and is not at all related to the properties of the trigger itself. In fact, the longest delays are due to simultaneous downlink sessions.

In attachment you can find the table used to produce the histogram. The columns are GRB name and corresponding delay time in seconds from the trigger. You can use it to generate the delays for the tests on the 1000-GRB sample. For the final sample, I can include this figure in the header of the fits files.

Best,

Maria Grazia

(atached file Swift_delay_times.txt)

Il 02/03/20 17:40, STOLARCZYK Thierry ha scritto:
> Dear Maria Grazia,
> Thanks a lot for this information.
> Do I understand that there is no correlation between the delay and the luminosity of the GRB ?
> Therefore the plot you show is also representative of the BRIGHT GRBs that are in the 1000-GRB sample ?
> If this is the case then can you give the raw data (list of numbers) or digitize the histogram so that we can generate the delays on the fly ?
> Thank you,
> Best regards,
> Thierry S.

> -----Message d'origine-----
> De : Maria Grazia Bernardini <maria.bernardini@inaf.it> Envoyé : 
> vendredi 28 février 2020 17:49 À : Giancarlo Ghirlanda 
> <giancarlo.ghirlanda@inaf.it>; Antonio Stamerra 
> <antonio.stamerra@inaf.it>; Francesco Longo 
> <francesco.longo@ts.infn.it>; Alessandro Carosi 
> <carosi@lapp.in2p3.fr>; SCHUSSLER Fabian <fabian.schussler@cea.fr>; 
> STOLARCZYK Thierry <thierry.stolarczyk@cea.fr>; Lara Nava 
> <lara.nava@inaf.it>; zeljka.bosnjak@fer.hr; 
> thomas.gasparetto@ts.infn.it; Iftach Sadeh <iftach.sadeh@desy.de> Cc : 
> Paolo D'Avanzo <paolo.davanzo@inaf.it>; Andrea Melandri 
> <andrea.melandri@inaf.it> Objet : delay time in Swift alerts
>
> Dear all,
>
> one of the issues raised during the last F2F meeting was to account for the typical delay in the delivery of the alerts for Swift. This delay is reported to be of the order of 20 s (see e.g. here:
> https://swift.gsfc.nasa.gov/analysis/suppl_uguide/GCN_mssg_info.html).
>
> To have a more realistic idea of this delay to include it in our 
> simulations, I checked the GRBs triggered by Swift from January 2017
> (238 events) and I have found that:
> - the minimum delay is 12 s;
> - the mean value is 77 s, and the median is 34 s;
> - 65% of GRBs have delays shorter than 52 s, 90% shorter than 130 s.
>
> You can see in attachment the distribution of these delays. I checked manually the bunch of cases where delays are very long (a few minutes), and usually this is due to the trigger occurring during a Malindi downlink session.
>
> Cheers,
>
> Maria Grazia

(attached : Swift_latency.eps)
