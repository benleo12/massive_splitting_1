(* Numerical evaluator for HyperFORM L() iterated integrals.
   Convention: L(a_1, ..., a_{n-1}, x) = G(a_1, ..., a_{n-1}; x).
   G(; x) = 1.
   G(a_1, ..., a_n; x) = \int_0^x dt/(t-a_1) G(a_2, ..., a_n; t).

   Shuffle regularization: leading 0s are peeled off via
     G(0^m b_1...b_n; x) = (log^m(x)/m!) G(b_1...b_n; x)
                          - sum over shuffles inserting 0s into the b-sequence.

   Caching: same L() value at the same x is computed once. *)

ClearAll[Lnum, shuffles];

shuffles[{}, b_List] := {b};
shuffles[a_List, {}] := {a};
shuffles[a_List, b_List] := Join[
  Map[Prepend[#, First[a]] &, shuffles[Rest[a], b]],
  Map[Prepend[#, First[b]] &, shuffles[a, Rest[b]]]
];

Lnum[{}, x_?NumericQ] := 1;

Lnum[lets_List, x_?NumericQ] := Lnum[lets, x] = Module[
  {m, rest, allShuffles, leadingZerosList, validShuffles},
  m = LengthWhile[lets, # === 0 &];
  rest = Drop[lets, m];
  Which[
    m == 0,
       NIntegrate[1/(s - First[lets]) Lnum[Rest[lets], s], {s, 0, x},
         AccuracyGoal -> 12, PrecisionGoal -> 8, WorkingPrecision -> 20,
         Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule"}],
    rest === {},
       Log[x]^m/m!,
    True,
       leadingZerosList = ConstantArray[0, m];
       allShuffles = shuffles[leadingZerosList, rest];
       validShuffles = DeleteCases[allShuffles, leadingZerosList ~Join~ rest];
       Log[x]^m/m! Lnum[rest, x] - Sum[Lnum[s, x], {s, validShuffles}]
  ]
];
