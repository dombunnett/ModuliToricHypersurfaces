# Bunnett: Computes the Cox ring and Automorphism group of a toric orbifold.
#
# First will compute the compute the class group
#

use strict;
use warnings;
use application "polytope";
use application "fulton";

##########################################
#FUNCTION TO MAKE MULTISETS
##########################################
sub make_multisets {
   my($m, $n) = @_;
   if($n == 1){
      return new Set<Vector<Integer>>([new Vector<Integer>([$m])]);
   }
   my $result = new Set<Vector<Integer>>();
   for(my $i=0; $i<$m+1; $i++){
      my $recurse = make_multisets($m-$i, $n-1);
      foreach my $v (@$recurse){
         $result += new Vector<Integer>([$i, @$v]);
      }
   }
   my $result2 = new Matrix<Integer>($result);
   return $result2;
}
#########################################

my $p = weighted_projective_space(new Vector<Int>(1,1,3));
my $n = $p->N_RAYS;
my $A =  new Matrix<Rational>(primitive($p->RAYS));

print $A,"\n";

my $partition = new Map<Int, Set<Int>>();
my $rays = new Set<Int>(sequence(0,$n));

# Next we loop over all the rays and check if they are equivalent, if
# they are, we remove them from the list and through it in a set.
# This is done via the exact sequence:
# 0 -> M -A-> ZZ^{#rays} -> Cl(X) -> 0.
# We check if two rays differ by an element in the image of A.

while($rays->size > 0){   
      my $first = $rays->front;
      my @class = grep{
            my $b = unit_vector($n,$_) - unit_vector($n,$first);
            my $test = new Polytope(EQUATIONS=>(-$b|$A), INEQUALITIES=>[unit_vector($n,0)]);
            $test->FEASIBLE;
      } @$rays;

      # print join(",",@class),"\n";
      # This prints all the members of the ray set which have the
      # same class as the one tested (in a reasonable way, hence
      # the join).
      $partition->{$first} = new Set<Int>(\@class);
      $rays -= $partition->{$first};
}
print $partition,"\n";

# The reductive part of the (an) automorphism group is given by
# the product of a general linear group for each subset in the
# partition. Hence we square the sizes and add for the dimension.

my $rays2 = new Set<Int>(sequence(0,$n));
my $red_dim = 0;
while($rays2->size >0){
      $red_dim = $red_dim + ($partition->{$rays2->front}->size)**2;
      $rays2 -=$rays2->front;
}

print "dimension of reductive part: ", $red_dim, "\n";


# Now to get the unipotent part we need to check for more 
# complicated relations between the unit vectors.
#
# We will check for example for [e_i + e_j] = [e_k].
# This will correspnd to monomials such as deg(x*y) = deg(z).

my $unip = new Map<Int, Int>();
for(my $i =0; $i<$n; $i++){
      $unip->{$i} = 0;
   for(my $m = 2; $m<10; $m++){
      # This 'if' below checks if the variable is of another class.
      # For example in PP(1,1,2) with coordinates x,y,z, we do not
      # need to check both x and y.
      if($partition->{$i}->size>0){
         my $rays = new Set<Int>(sequence(0,$n));
         $rays -= $partition->{$i};
         my $N = $rays->size;
         
         if($N==1){
            for(my $r=1; $r<10; $r++){
               my $j = @$rays[0];
               my $b = unit_vector($n,$i) - $r*unit_vector($n,$j);
               my $test = new Polytope(EQUATIONS=>(-$b|$A), INEQUALITIES=>[unit_vector($n,0)]);
               if($test->FEASIBLE){
                  $unip->{$i} += 1;
               }
            }
         }
         if($N > 1){
            
            my $M = make_multisets($m,$N);

            # my $B = unit_matrix($N);
            # my $P = new Polytope(INEQUALITIES=>zero_vector|$B,EQUATIONS=>-$m|ones_vector($N));
            # my $M = $P->LATTICE_POINTS;
            #
            # $P is the correct simplex to give us the multisets.
            # 
            # This first for loop loops over the multisets. Thus
            # if we find a FEASIBLE vector, it should go in a set.
            for(my $l=0; $l < $M->rows; $l++){

               # $v will be our non-variable monomial.
               my $v = unit_vector($n,1)-unit_vector($n,1);
               for(my $k=0; $k<$N;$k++){
                  my $var = @$rays[$k];
                  $v += $M->[$l][$k]*unit_vector($n,$var);
               }
               # now we can check whether or not $v is the
               # same as any of the other variables.
               my $b = $v - unit_vector($n,$i);
               my $test = new Polytope(EQUATIONS=>(-$b|$A), INEQUALITIES=>[unit_vector($n,0)]);
               if($test->FEASIBLE){
                  # HERE I JUST COUNT THE NUMBER OF MONOMIALS
                  # IT WOULD BE GOOD TO ACTUALLY RECORD THE 
                  # MONOMIAL ITSELF AND NOT JUST THE NUMBER OF
                  # THEM.
                  $unip->{$i} += 1;
               }
            }
         }
      }
   }
}
print $unip, "\n";

# Now we actually compute the dimension of the unipotent part
# using the numbers above by a result of Cox.

my $unip_dim = 0;
for(my $i=0; $i<$n; $i++){
   $unip_dim += ($unip->{$i})*($partition->{$i}->size);
   print "ray: ", $i, ", unipotent roots: ", $unip->{$i}, ", number of equiv. vars: ", $partition->{$i}->size, "\n";
}

print "dimension of the unipotent part: ", $unip_dim,"\n";

# Finally, the dimension is just the sum of both parts!

declare $dimension = $red_dim + $unip_dim;
print "Dimension of the graded automoprhism group is: ", $dimension;
