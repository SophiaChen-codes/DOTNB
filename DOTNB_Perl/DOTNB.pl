use strict;
use warnings;
use constant pi => 4*atan2(1, 1);

sub Gamma { #gamma function using lanczos method
    my $z = shift;
    use constant g => 9;
    $z <  .5 ?  pi / sin(pi * $z) / Gamma(1 - $z) :
    sqrt(2* pi) * ($z + g - .5)**($z - .5) * exp(-($z + g - .5)) *
    do {
        my @coeff = qw{
            1.000000000000000174663
            5716.400188274341379136
            -14815.30426768413909044
            14291.49277657478554025
            -6348.160217641458813289
            1301.608286058321874105
            -108.1767053514369634679
            2.605696505611755827729
           -0.7423452510201416151527e-2
           0.5384136432509564062961e-7
           -0.4023533141268236372067e-8
        };
        my ($sum, $i) = (shift(@coeff), 0);
        $sum += $_ / ($z + $i++) for @coeff;
        $sum;
    }
} 

sub isinf { $_[0]==9**9**9 || $_[0]==-9**9**9 }

sub Hypergeom {
## hypergeometric function
    my $a=shift;
    my $b=shift;
    my $c=shift;
    my $x=shift;
    my $TOLERANCE=1.0e-10;
    my $term=$a * $b * $x / $c; 
    my $value =1.0 + $term;
    my $n=1;
    while( abs($term) > $TOLERANCE){
        $a++; $b++; $c++; $n++;
        $term*=$a*$b*$x / $c /$n;
        $value+= $term;
    }
    $value;
}

#printf "gamma(4)=%.12f\n", Gamma(4);
#printf "Hypergeom(2,3,4,0.5)=%.12f\n", Hypergeom(2,3,4,0.5);

sub DOTNB_pdf{
  # calculate the p-value of two different NB distributions (r1,p1) and (r2,p2)
  # dn is the observed number of difference
    my $r1=shift;   my $p1=shift;
    my $r2=shift;   my $p2=shift;
    my $dn=shift;
    my $q1 = 1-$p1; my $q2 = 1-$p2;
	my $prob;
    if($dn>0){
        $prob = ($p1**$r1)*($p2**$r2)*((Gamma($r1+$dn)/Gamma($r1))/Gamma(1+$dn))*($q1**$dn)*Hypergeom($r1+$dn,$r2,$dn+1,$q1*$q2);
    }else{
        $dn=-$dn;
        $prob = ($p1**$r1)*($p2**$r2)*((Gamma($r2+$dn)/Gamma($r2))/Gamma(1+$dn))*($q2**$dn)*Hypergeom($r2+$dn,$r1,$dn+1,$q1*$q2);
    }
    $prob;
}

sub DOTNB_mean{
  # calculate the mean of the difference of two NB distributions (r1,p1) and (r2,p2)
	my $r1=shift;   my $p1=shift;
    my $r2=shift;   my $p2=shift;
	my $q1=1-$p1;   my $q2=1-$p2;
	my $dmean=(($r1*$q1)/$p1) - ($r2*$q2/$p2);
	$dmean;
}

sub DOTNB_var{
  # calculate the variance of the difference of two NB distributions (r1,p1) and (r2,p2)
    my $r1=shift;   my $p1=shift;
    my $r2=shift;   my $p2=shift;
    my $q1 = 1-$p1;
    my $q2 = 1-$p2;
    my $dvar = ($r1*$q1)/($p1^2) + ($r2*$q2)/($p2^2);
    $dvar;
}

sub DOTNB_cdf{
  # calculate the cumulative probability of dn that is the difference of two NB distributions (r1,p1) and (r2,p2)
  # dn is the observed number of difference
  # sum pvalue for <= dn, 
	my $r1=shift;   my $p1=shift;
    my $r2=shift;   my $p2=shift;
    my $dn=shift; 
	my $cvalue=0;
	my $pdfv=DOTNB_pdf($r1,$p1,$r2,$p2,$dn);
	if(isinf($pdfv)) {
		if ($dn > DOTNB_mean($r1,$p1,$r2,$p2)){
			$cvalue=1;
		}
	}else{
		while(not isinf($pdfv)){
			$cvalue=$cvalue+$pdfv;
			$dn=$dn-1;
			$pdfv=DOTNB_pdf($r1,$p1,$r2,$p2,$dn);
		}
	}
	$cvalue;
}

#test
my $r1=2; my $p1=0.58; my $r2=3; my $p2=0.3,my $k=100;
print "DOTNB_pdf(r1=$r1, p1=$p1, r2=$r2, p2=$p2, dn=$k)=", DOTNB_pdf($r1,$p1,$r2,$p2,$k), "\n";
print "DOTNB_cdf(r1=$r1, p1=$p1, r2=$r2, p2=$p2, dn=$k)=", DOTNB_cdf($r1,$p1,$r2,$p2,$k), "\n";


