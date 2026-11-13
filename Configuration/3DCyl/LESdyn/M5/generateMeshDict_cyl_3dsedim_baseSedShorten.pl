#!/usr/bin/perl
# Generate the cylinder blockMeshDict file for the von Karman LES case
#              16          17         18          19
#  y ^        0 - - - - - 1- - - - - - 2 - - - - - 3
#    |        |           |            |           |  
#    --> x    |     0     |    1       |   2       |     ny1
#             |20         |21          |22         |23
#             4 - - - - - 5 - - - - - -6 - - - - - 7
#
#             |           |32/36 33/37 |           |  
#
#             |    3      |34/38 35/39 |   4       |     ny2
#
#             |24         |25          |26         |27
#             8 - - - - - 9 - - - - - 10 - - - - -11
#             |           |            |           |  
#             |    5      |   6        |   7       |     ny1
#             |28         |29          |30         |31
#             12 - - - - 13 - - - - - 14 - - - - - 15
#                                                   
#                 nx1       nx2=ny2       nx3
use strict;

if ($#ARGV != 8)
{
    usage();
    exit;
}

my @input = @ARGV;
my $D = $ARGV[0];
my $H0 = $ARGV[1];
my $H1 = $ARGV[2];
my $H1b = $ARGV[3];
my $H2 = $ARGV[4];
my $H3 = $ARGV[5]; #for the bottom block
my $L1 = $ARGV[6];
my $L2 = $ARGV[7];
my $W = $ARGV[8];

my $w = 0.5*$W;
my $r = 0.5*$D;
my $r0 = 0.5*$D;
my $rinner1 = $r0/4.;
my $rinner = $rinner1/sqrt(2.);
my $xe = $L2-$L1;
my $rext = $r*sqrt(2.);


# aspect ratio for largest/smallest cells for the different blocks M0123456
my $g1 = 1.;
my $g2 = 2.3632098693316084;#1.1511012969725847 1.5456502002509844 1.6508780037043267 2.037381869748463 2.3632098693316084 2.6590599489095053 2.9109690503215853
my $g3 = 3.9665682792349184;#3.630890385668341 4.075542450220805 3.856223782579019 4.298316982343356 4.103832726447115 3.9665682792349184 4.263778908799975
my $gC = 2.522311151573334;#2.7416556078516723 2.4448365874972424 2.5854939547329243 2.3292823634703046 2.43820920811315 2.522311151573334 2.3447199897377486
my $gz = 1.;
my $gZ = 112.1832662895043;#41.432951063305104 62.33183979817063 62.33183979817063 112.1832662895043 112.1832662895043 112.1832662895043 190.15441154173178
my $gZ1 = 0.0775727693174146;#0.5830419509076277 0.21697282450595096  0.21697282450595096 0.0775727693174146 0.0775727693174146 0.0775727693174146 0.04204588825765332
my $gZ1b = 0.9351357050411135;#0.2652865369156223 0.5419339335730116 0.5419339335730116 0.9351357050411135 0.9351357050411135 0.9351357050411135 0.954366129448506

my $nx1 = 23;#7 10 13 17 20 23 27
my $nx2 = 70;#22 30 40 50 60 70 80
my $nx3 = 47;#15 20 27 33 40 47 53

my $ny1 = 12;#6 7 9 10 11 12 13
my $ny2 = $nx2;

my $nz = 12;#11 11 11 12 12 12 14
my $nzb = 58;#42 46 46 58 58 58 76
my $nz1 = 150;#176 148 148 150 150 150 194
my $nz1b = 42;#52 44 44 42 42 42 46

my $g1i = 1.0/$g1;
my $g2i = 1.0/$g2;
my $g3i = 1.0/$g3;
my $gCi = 1.0/$gC;
my $gZi = 1.0/$gZ;
my $gZi1 = $gZ1;
my $gZi1b = $gZ1b;

printHeader("2.0");
printConvertToMeters("1.0");
printVertices();
printBlocks();
printEdges();
printPatches();
printMerge();

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

sub printHeader
{
    print "FoamFile\n";
    print "{\n";
    print "    version\t $_[0];\n";
    print "    format\t ascii;\n";
    print "    class\t dictionary;\n";
    print "    object\t blockMeshDict;\n";
    print "}\n\n";

    print "/*\n";
    print "Automatically generated with $0 @input\n";
    print "*/\n\n";
}

sub printConvertToMeters
{
    print "convertToMeters\t $_[0];\n\n";
}

sub printVertices
{
    print "vertices\n";
    print "(\n";
    print "    ( -$L1  $w $H1 ) //0\n";
    print "    ( -$r   $w $H1 ) //1\n";
    print "    (  $r   $w $H1 ) //2\n";
    print "    (  $xe  $w $H1 ) //3\n";

    print "    ( -$L1  $r $H1 ) //4\n";
    print "    ( -$r   $r $H1 ) //5\n";
    print "    (  $r   $r $H1 ) //6\n";
    print "    (  $xe  $r $H1 ) //7\n";

    print "    ( -$L1  -$r $H1 ) //8\n";
    print "    ( -$r   -$r $H1 ) //9\n";
    print "    (  $r   -$r $H1 ) //10\n";
    print "    (  $xe  -$r $H1 ) //11\n";

    print "    ( -$L1  -$w $H1 ) //12\n";
    print "    ( -$r   -$w $H1 ) //13\n";
    print "    (  $r   -$w $H1 ) //14\n";
    print "    (  $xe  -$w $H1 ) //15\n";

    print "    ( -$L1  $w $H1b ) //16\n";
    print "    ( -$r   $w $H1b ) //17\n";
    print "    (  $r   $w $H1b ) //18\n";
    print "    (  $xe  $w $H1b ) //19\n";

    print "    ( -$L1  $r $H1b ) //20\n";
    print "    ( -$r   $r $H1b ) //21\n";
    print "    (  $r   $r $H1b ) //22\n";
    print "    (  $xe  $r $H1b ) //23\n";

    print "    ( -$L1  -$r $H1b ) //24\n";
    print "    ( -$r   -$r $H1b ) //25\n";
    print "    (  $r   -$r $H1b ) //26\n";
    print "    (  $xe  -$r $H1b ) //27\n";

    print "    ( -$L1  -$w $H1b ) //28\n";
    print "    ( -$r   -$w $H1b ) //29\n";
    print "    (  $r   -$w $H1b ) //30\n";
    print "    (  $xe  -$w $H1b ) //31\n";
    
    print "    ( -$L1  $w $H2 )   //32\n";
    print "    ( -$r   $w $H2 )   //33\n";
    print "    (  $r   $w $H2 )   //34\n";
    print "    (  $xe  $w $H2 )   //35\n";

    print "    ( -$L1  $r $H2 )   //36\n";
    print "    ( -$r   $r $H2 )   //37\n";
    print "    (  $r   $r $H2 )   //38\n";
    print "    (  $xe  $r $H2 )   //39\n";

    print "    ( -$L1  -$r $H2 )  //40\n";
    print "    ( -$r   -$r $H2 )  //41\n";
    print "    (  $r   -$r $H2 )  //42\n";
    print "    (  $xe  -$r $H2 )  //43\n";

    print "    ( -$L1  -$w $H2 )  //44\n";
    print "    ( -$r   -$w $H2 )  //45\n";
    print "    (  $r   -$w $H2 )  //46\n";
    print "    (  $xe  -$w $H2 )  //47\n";
    
# vertices on the cylinder
    print "    ( -$rinner   $rinner $H1 )  //48\n";
    print "    (  $rinner   $rinner $H1 )  //49\n";
    print "    ( -$rinner   -$rinner $H1 ) //50\n";
    print "    (  $rinner   -$rinner $H1 ) //51\n";

    print "    ( -$rinner   $rinner $H1b )  //52\n";
    print "    (  $rinner   $rinner $H1b )  //53\n";
    print "    ( -$rinner   -$rinner $H1b ) //54\n";
    print "    (  $rinner   -$rinner $H1b ) //55\n";
    
    print "    ( -$rinner   $rinner $H2 )  //56\n";
    print "    (  $rinner   $rinner $H2 )  //57\n";
    print "    ( -$rinner   -$rinner $H2 ) //58\n";
    print "    (  $rinner   -$rinner $H2 ) //59\n";

# Lower blocks
    print "    ( -$L1  $w $H0 ) //60\n";
    print "    ( -$r   $w $H0 ) //61\n";
    print "    (  $r   $w $H0 ) //62\n";
    print "    (  $xe  $w $H0 ) //63\n";

    print "    ( -$L1  $r $H0 ) //64\n";
    print "    ( -$r   $r $H0 ) //65\n";
    print "    (  $r   $r $H0 ) //66\n";
    print "    (  $xe  $r $H0 ) //67\n";

    print "    ( -$L1  -$r $H0 ) //68\n";
    print "    ( -$r   -$r $H0 ) //69\n";
    print "    (  $r   -$r $H0 ) //70\n";
    print "    (  $xe  -$r $H0 ) //71\n";

    print "    ( -$L1  -$w $H0 ) //72\n";
    print "    ( -$r   -$w $H0 ) //73\n";
    print "    (  $r   -$w $H0 ) //74\n";
    print "    (  $xe  -$w $H0 ) //75\n";

    print "    ( -$rinner   $rinner $H0 ) //76\n";
    print "    (  $rinner   $rinner $H0 ) //77\n";
    print "    ( -$rinner   -$rinner $H0 ) //78\n";
    print "    (  $rinner   -$rinner $H0 ) //79\n";


# Scour pit blocks

    print "    ( -$r   $r $H3 ) //80\n";
    print "    (  $r   $r $H3 ) //81\n";

    print "    ( -$r   -$r $H3 ) //82\n";
    print "    (  $r   -$r $H3 ) //83\n";

    print "    ( -$rinner   $rinner $H3 ) //84\n";
    print "    (  $rinner   $rinner $H3 ) //85\n";
    print "    ( -$rinner   -$rinner $H3 ) //86\n";
    print "    (  $rinner   -$rinner $H3 ) //87\n";

    print ");\n\n";
}

sub printBlocks
{
    print "blocks\n";
    print "(\n";
    print "    hex (4 5 1 0 20 21 17 16) ($nx1 $ny1 $nzb) simpleGrading ($g1i $g2 $gZ) //0\n";
    print "    hex (5 6 2 1 21 22 18 17) ($nx2 $ny1 $nzb) simpleGrading (1 $g2 $gZ) //1\n";
    print "    hex (6 7 3 2 22 23 19 18) ($nx3 $ny1 $nzb) simpleGrading ($g3 $g2 $gZ) //2\n";

    print "    hex (8 9 5 4 24 25 21 20) ($nx1 $ny2 $nzb) simpleGrading ($g1i 1 $gZ) //3\n";
    print "    hex (10 11 7 6 26 27 23 22) ($nx3 $ny2 $nzb) simpleGrading ($g3 1 $gZ)//4\n";

    print "    hex (12 13 9 8 28 29 25 24) ($nx1 $ny1 $nzb) simpleGrading ($g1i $g2i $gZ)//5\n";
    print "    hex (13 14 10 9 29 30 26 25) ($nx2 $ny1 $nzb) simpleGrading (1 $g2i $gZ)//6\n";
    print "    hex (14 15 11 10 30 31 27 26) ($nx3 $ny1 $nzb) simpleGrading ($g3 $g2i $gZ)//7\n";

    print "    hex (48 49 6 5 52 53 22 21) ($nx2 $ny2 $nzb) simpleGrading (1 $gC $gZ)//8\n";
    print "    hex (51 10 6 49 55 26 22 53) ($nx2 $ny2 $nzb) simpleGrading ($gC 1 $gZ)//9\n";
    print "    hex (9 10 51 50 25 26 55 54) ($nx2 $ny2 $nzb) simpleGrading (1 $gCi $gZ)//10\n";
    print "    hex (9 50 48 5 25 54 52 21) ($nx2 $ny2 $nzb) simpleGrading ($gCi 1 $gZ)//11\n";
    
    print "    hex (20 21 17 16 36 37 33 32) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2 $gz) //12\n";
    print "    hex (21 22 18 17 37 38 34 33) ($nx2 $ny1 $nz) simpleGrading (1 $g2 $gz) //13\n";
    print "    hex (22 23 19 18 38 39 35 34) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2 $gz) //14\n";

    print "    hex (24 25 21 20 40 41 37 36) ($nx1 $ny2 $nz) simpleGrading ($g1i 1 $gz) //15\n";
    print "    hex (26 27 23 22 42 43 39 38) ($nx3 $ny2 $nz) simpleGrading ($g3 1 $gz)//16\n";

    print "    hex (28 29 25 24 44 45 41 40) ($nx1 $ny1 $nz) simpleGrading ($g1i $g2i $gz)//17\n";
    print "    hex (29 30 26 25 45 46 42 41) ($nx2 $ny1 $nz) simpleGrading (1 $g2i $gz)//18\n";
    print "    hex (30 31 27 26 46 47 43 42) ($nx3 $ny1 $nz) simpleGrading ($g3 $g2i $gz)//19\n";

    print "    hex (52 53 22 21 56 57 38 37) ($nx2 $ny2 $nz) simpleGrading (1 $gC $gz)//20\n";
    print "    hex (55 26 22 53 59 42 38 57) ($nx2 $ny2 $nz) simpleGrading ($gC 1 $gz)//21\n";
    print "    hex (25 26 55 54 41 42 59 58) ($nx2 $ny2 $nz) simpleGrading (1 $gCi $gz)//22\n";
    print "    hex (25 54 52 21 41 58 56 37) ($nx2 $ny2 $nz) simpleGrading ($gCi 1 $gz)//23\n";

    print "    hex (64 65 61 60 4 5 1 0) ($nx1 $ny1 $nz1) simpleGrading ($g1i $g2 $gZ1) //24\n";
    print "    hex (65 66 62 61 5 6 2 1) ($nx2 $ny1 $nz1) simpleGrading (1 $g2 $gZ1) //25\n";
    print "    hex (66 67 63 62 6 7 3 2) ($nx3 $ny1 $nz1) simpleGrading ($g3 $g2 $gZ1) //26\n";

    print "    hex (68 69 65 64 8 9 5 4) ($nx1 $ny2 $nz1) simpleGrading ($g1i 1 $gZ1) //27\n";
    print "    hex (70 71 67 66 10 11 7 6) ($nx3 $ny2 $nz1) simpleGrading ($g3 1 $gZ1)//28\n";

    print "    hex (72 73 69 68 12 13 9 8) ($nx1 $ny1 $nz1) simpleGrading ($g1i $g2i $gZ1)//29\n";
    print "    hex (73 74 70 69 13 14 10 9) ($nx2 $ny1 $nz1) simpleGrading (1 $g2i $gZ1)//30\n";
    print "    hex (74 75 71 70 14 15 11 10) ($nx3 $ny1 $nz1) simpleGrading ($g3 $g2i $gZ1)//31\n";

    print "    hex (76 77 66 65 48 49 6 5) ($nx2 $ny2 $nz1) simpleGrading (1 $gC $gZ1)//32\n";
    print "    hex (79 70 66 77 51 10 6 49) ($nx2 $ny2 $nz1) simpleGrading ($gC 1 $gZ1)//33\n";
    print "    hex (69 70 79 78 9 10 51 50) ($nx2 $ny2 $nz1) simpleGrading (1 $gCi $gZ1)//34\n";
    print "    hex (69 78 76 65 9 50 48 5) ($nx2 $ny2 $nz1) simpleGrading ($gCi 1 $gZ1)//35\n";

    print "    hex (84 85 81 80 76 77 66 65) ($nx2 $ny2 $nz1b) simpleGrading (1 $gC $gZ1b)//36\n";
    print "    hex (87 83 81 85 79 70 66 77) ($nx2 $ny2 $nz1b) simpleGrading ($gC 1 $gZ1b)//37\n";
    print "    hex (82 83 87 86 69 70 79 78) ($nx2 $ny2 $nz1b) simpleGrading (1 $gCi $gZ1b)//38\n";
    print "    hex (82 86 84 80 69 78 76 65) ($nx2 $ny2 $nz1b) simpleGrading ($gCi 1 $gZ1b)//39\n";


    print ");\n\n";
}

sub printEdges
{
    print "edges\n";
    print "(\n";
    print "    arc 48 49 ( 0 $rinner1 $H1 )\n";
    print "    arc 49 51 ( $rinner1 0 $H1 )\n";
    print "    arc 50 51 ( 0 -$rinner1 $H1 )\n";
    print "    arc 50 48 ( -$rinner1 0 $H1 )\n";
    
    print "    arc 52 53 ( 0 $rinner1 $H1b )\n";
    print "    arc 55 53 ( $rinner1 0 $H1b )\n";
    print "    arc 54 55 ( 0 -$rinner1 $H1b )\n";
    print "    arc 54 52 ( -$rinner1 0 $H1b )\n";

    print "    arc 56 57 ( 0 $rinner1 $H2 )\n";
    print "    arc 59 57 ( $rinner1 0 $H2 )\n";
    print "    arc 58 59 ( 0 -$rinner1 $H2 )\n";
    print "    arc 58 56 ( -$rinner1 0 $H2 )\n";

    print "    arc 76 77 ( 0 $rinner1 $H0 )\n";
    print "    arc 79 77 ( $rinner1 0 $H0 )\n";
    print "    arc 78 79 ( 0 -$rinner1 $H0 )\n";
    print "    arc 78 76 ( -$rinner1 0 $H0 )\n";
    
    print "    arc 84 85 ( 0 $rinner1 $H3 )\n";
    print "    arc 87 85 ( $rinner1 0 $H3 )\n";
    print "    arc 86 87 ( 0 -$rinner1 $H3 )\n";
    print "    arc 86 84 ( -$rinner1 0 $H3 )\n";

    print "    arc 5 6 ( 0 $rext $H1 )\n";
    print "    arc 6 10 ( $rext 0 $H1 )\n";
    print "    arc 9 10 ( 0 -$rext $H1 )\n";
    print "    arc 9 5 ( -$rext 0 $H1 )\n";
    
    print "    arc 21 22 ( 0 $rext $H1b )\n";
    print "    arc 22 26 ( $rext 0 $H1b )\n";
    print "    arc 25 26 ( 0 -$rext $H1b )\n";
    print "    arc 25 21 ( -$rext 0 $H1b )\n";

    print "    arc 37 38 ( 0 $rext $H2 )\n";
    print "    arc 38 42 ( $rext 0 $H2 )\n";
    print "    arc 41 42 ( 0 -$rext $H2 )\n";
    print "    arc 41 37 ( -$rext 0 $H2 )\n";

    print "    arc 65 66 ( 0 $rext $H0 )\n";
    print "    arc 66 70 ( $rext 0 $H0 )\n";
    print "    arc 69 70 ( 0 -$rext $H0 )\n";
    print "    arc 69 65 ( -$rext 0 $H0 )\n";

    print "    arc 80 81 ( 0 $rext $H3 )\n";
    print "    arc 83 81 ( $rext 0 $H3 )\n";
    print "    arc 82 83 ( 0 -$rext $H3 )\n";
    print "    arc 82 80 ( -$rext 0 $H3 )\n";


    print ");\n\n";
}

sub printPatches
{
    print "boundary\n";
    print "(\n";
    print "    cylinder\n ";
    print "    { \n";
    print " 	type wall;\n";
    print "        faces \n";
    print "        (\n";
    print "        (48 52 53 49)\n";
    print "        (49 53 55 51)\n";
    print "        (51 55 54 50)\n";
    print "        (50 54 52 48)\n";
    print "        (52 56 57 53)\n";
    print "        (53 57 59 55)\n";
    print "        (55 59 58 54)\n";
    print "        (54 58 56 52)\n";
    print "        (76 48 49 77)\n";
    print "        (77 49 51 79)\n";
    print "        (79 51 50 78)\n";
    print "        (78 50 48 76)\n";
    print "        (84 76 77 85)\n";
    print "        (85 77 79 87)\n";
    print "        (87 79 78 86)\n";
    print "        (86 78 76 84)\n";
    print "        );\n";
    print "    }\n\n";

    print "    inlet\n";
    print "    {\n";
    print "	type patch;\n";
    print "        faces\n";
    print "        (\n";
    print "        (4 20 16 0)\n";
    print "        (8 24 20 4)\n";
    print "        (12 28 24 8)\n";
    print "        (20 36 32 16)\n";
    print "        (24 40 36 20)\n";
    print "        (28 44 40 24)\n";
    print "        (64 4 0 60)\n";
    print "        (68 8 4 64)\n";
    print "        (72 12 8 68)\n";
    print "        );\n";
    print "    }\n\n";
       
    print "    outlet\n";
    print "    {\n";
    print "	type patch;\n";
    print "        faces\n";
    print "        (\n";
    print "        (3 19 23 7)\n";
    print "        (7 23 27 11)\n";
    print "        (11 27 31 15)\n";
    print "        (19 35 39 23)\n";
    print "        (23 39 43 27)\n";
    print "        (27 43 47 31)\n";
    print "        );\n";
    print "    }\n\n";


    print "    outletb\n";
    print "    { \n";
    print "	type wall;\n";
    print "        faces \n";
    print "        (\n";
    print "        (63 3 7 67)\n";
    print "        (67 7 11 71)\n";
    print "        (71 11 15 75)\n";
    print "        );\n";
    print "    }\n\n";


    print "    pit\n";
    print "    { \n";
    print "	type wall;\n";
    print "        faces \n";
    print "        (\n";
    print "        (65 66 81 80)\n";
    print "        (82 83 70 69)\n";
    print "        (81 66 70 83)\n";
    print "        (82 69 65 80)\n";
    print "        );\n";
    print "    }\n\n";


    print "    lateral\n";
    print "    {\n";
    print "	type cyclic;\n";
    print "	neighbourPatch symplane;\n";
    print "	matchTolerance 0.0001;\n";
    print "        faces\n";
    print "        (\n";
    print "        (16 17 1 0)\n";
    print "        (17 18 2 1)\n";
    print "        (18 19 3 2)\n";
    print "        (32 33 17 16)\n";
    print "        (33 34 18 17)\n";
    print "        (34 35 19 18)\n";
    print "        (0 1 61 60)\n";
    print "        (1 2 62 61)\n";
    print "        (2 3 63 62)\n";
    print "        );\n";
    print "    }\n\n";

    print "    symplane\n";
    print "    {\n";
    print "	type cyclic;\n";
    print "     neighbourPatch lateral;\n";
    print "	matchTolerance 0.0001;\n";
    print "        faces\n";
    print "        (\n";
    print "        (12 13 29 28)\n";
    print "        (13 14 30 29)\n";
    print "        (14 15 31 30)\n";
    print "        (28 29 45 44)\n";
    print "        (29 30 46 45)\n";
    print "        (30 31 47 46)\n";
    print "        (72 73 13 12)\n";
    print "        (73 74 14 13)\n";
    print "        (74 75 15 14)\n";
    print "        );\n";
    print "    }\n\n";

    print "    bottom\n";
    print "    { \n";
    print "	type wall;\n";
    print "        faces \n";
    print "        (\n";
    print "        (60 61 65 64)\n";
    print "        (61 62 66 65)\n";
    print "        (62 63 67 66)\n";
    print "        (64 65 69 68)\n";
    print "        (66 67 71 70)\n";
    print "        (68 69 73 72)\n";
    print "        (69 70 74 73)\n";
    print "        (70 71 75 74)\n";
    print "        (80 81 85 84)\n";
    print "        (81 83 87 85)\n";
    print "        (86 87 83 82)\n";
    print "        (80 84 86 82)\n";
    print "        );\n";
    print "    }\n\n";

    print "    surface\n";
    print "    {\n";
    print "	type patch;\n";
    print "        faces\n";
    print "        (\n";
    print "        (36 37 33 32)\n";
    print "        (37 38 34 33)\n";
    print "        (38 39 35 34)\n";
    print "        (40 41 37 36)\n";
    print "        (42 43 39 38)\n";
    print "        (44 45 41 40)\n";
    print "        (45 46 42 41)\n";
    print "        (46 47 43 42)\n";
    print "        (56 57 38 37)\n";
    print "        (59 42 38 57)\n";
    print "        (41 42 59 58)\n";
    print "        (41 58 56 37)\n";
    print "        );\n";
    print "    }\n\n";
    print ");\n\n";
}

sub printMerge
{
    print "mergPatchPairs\n";
    print "(\n";
    print ");\n";
}

sub usage
{
    print "usage:\n";
    print "\t $0 <D> <H0> <H1> <H1b> <H2> <H3> <L1> <L2> <W> <delta>\n";
    print "\t where (recommended value):\n";
    print "\t D   = Width of the refinment region around the cylinder. (4)\n";
    print "\t H0  = Elevation of the fixed sediment bottom. (-0.25)\n";
    print "\t H1  = Elevation of the sediment bed-fluid interface. (0.)\n";
    print "\t H1b = Elevation of the sediment bed-fluid interface. (1.)\n";    
    print "\t H2  = Elevation of the top boundary. (2.)\n";
    print "\t H3  = Elevation of the pit boundary. (-1.)\n";
    print "\t L1  = Length from inlet to cylinder center. (4.)\n";
    print "\t L2  = Length of computational domain. (12.)\n";
    print "\t W   = Width of the domain (8)\n";
    print "\t delta = cell size across the sides of the cylinder. (0.0023)\n";
}

