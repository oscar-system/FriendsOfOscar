using S_characters
using Oscar

# the smallest example of an S-character
# that is nonzero on all classes of prime power order elements
tbl = character_table("A8");
res = s_characters(t, rational = false, ppow_nonzero = true)
