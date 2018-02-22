Checklist for ...:

  -- Check everything is intialised properly - if in doubt initialise EVERYTHING to
  zero (in particular, check youve passed tensors the right number of values FFS)

  -- If the initial maximum residual changes each time you run the code then
  initial data isnt the same with each run.

  --Set so the initial data is forced to always be the same - are the initial
  maximum residuals the same with each run now? If not, there's a problem
  with the residual calculation

  -- If no newton iterations changes with each run, then there may be a problem
  in the Jacobian calculation. - comment out fill_in_contribution_to_jacobian to
  finite difference the jacobian. If it works now then there's a problem in the
  jacobian calculation 
