#!/usr/bin/perl -w


#Only calculate two highest allel genotypes
#2019年 12月 27日 星期五 15:00:03 CST
#
#新改动，每测到的地方改用-8作为没测到的标记。算LOD时，如果一个windou里有一半（暂定，根据剩下多少来算）都是-8，就把这个window从统计中舍弃

#Tue Dec 13 15:05:11 CST 2011


#eg of input
#chromosome07    2871346 G       G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G G - - - G - G - G G - G G G - G - G - G G G - -


#eg of output
#			ind	jap	niv	ruf	wild	ind_niv	jap_ruf	ind_wild	jap_wild
#chr		2312	0	0	0	0	0	0	0	0		0

$marker = 0;

$fname = $ARGV[0];
open OUT, ">$fname.pi" or die;

my %Abbrev = (
                'A' => [ 'A' ],
                'C' => [ 'C' ],
                'G' => [ 'G' ],
                'T' => [ 'T' ],
                'M' => [ 'A', 'C' ],
                'R' => [ 'A', 'G' ],
                'W' => [ 'A', 'T' ],
                'S' => [ 'C', 'G' ],
                'Y' => [ 'C', 'T' ],
                'K' => [ 'G', 'T' ],
                'V' => [ 'A', 'C', 'G' ],
                'H' => [ 'A', 'C', 'T' ],
                'D' => [ 'A', 'G', 'T' ],
                'B' => [ 'C', 'G', 'T' ],
                'X' => [ 'A', 'C', 'G', 'T' ],
                'N' => [ 'A', 'C', 'G', 'T' ]
);


#1-12    indica
#13-35   japonica
#36-45   nivara
#46-60   rufipogon

my @group_aus = qw/13 15 86 217 222 227 230 236 250 251 271 274 278 292 293 294 296 297 298 300 301 302 305 306 307 308 310 317 323 333 336 342 343 344 414 415 476 477 478 479 485 486 491 495 496 529/;
my @group_japonica = qw/3 4 5 10 12 14 16 17 18 23 26 28 29 32 34 35 48 52 54 56 57 63 64 67 70 71 74 75 79 82 83 84 85 87 89 93 101 103 106 111 116 119 120 123 129 132 133 136 137 143 145 148 149 151 152 153 170 171 176 178 180 184 185 186 187 188 194 195 204 205 206 207 208 209 210 212 213 214 215 218 219 220 221 223 224 226 229 240 243 247 249 254 255 258 259 261 262 263 264 265 266 267 269 270 284 290 304 311 312 315 318 319 320 321 322 324 325 329 332 341 345 362 364 375 394 407 408 409 425 426 448 449 450 452 453 456 459 472 474 480 482 483 484 487 489 499 501 502 504 518 522 523 524 525 526 528/;

#indica acutally
my @group_indica = qw/1 6 7 8 9 11 19 20 21 22 24 25 27 30 31 33 36 37 38 39 40 41 42 43 44 45 46 47 49 50 53 55 58 59 60 61 62 65 66 68 69 72 73 76 77 78 80 81 88 90 91 92 94 95 96 97 98 99 100 102 104 105 107 108 109 110 112 113 114 115 117 118 121 122 124 125 126 127 128 130 131 134 135 138 139 140 141 142 144 146 147 150 154 155 156 157 158 159 160 161 162 163 164 165 166 167 168 169 172 173 174 175 177 179 181 182 183 189 190 191 192 193 196 197 198 199 200 201 202 203 211 225 228 231 232 233 234 235 237 238 241 242 245 246 248 253 256 260 268 272 273 275 276 277 280 282 285 286 287 288 289 291 295 303 313 314 316 326 328 334 335 339 340 346 347 348 349 350 351 352 353 354 355 356 357 358 359 361 363 365 366 367 368 369 370 371 372 373 374 376 378 379 380 381 382 383 384 385 386 387 388 389 390 391 395 396 397 398 399 400 401 402 403 404 405 406 410 411 412 413 416 417 418 421 422 423 424 427 428 429 430 431 432 433 434 435 436 437 438 439 440 441 442 443 444 445 446 447 458 461 462 463 464 465 466 467 468 469 470 471 473 475 481 488 490 492 494 497 498 500 506 507 508 509 510 511 512 513 514 515 516 517 519 520 521/;

#print (int(@group_aus), "\t", int(@group_japonica), "\t", int(@group_indica), "\n");
#exit;

my @group_aus_indica = @group_aus;
push(@group_aus_indica, @group_indica);
my @group_jap_indica = @group_japonica;
push(@group_jap_indica, @group_indica);

$minimum_accession = 5;


my $print_cult = "";
my $print_wild = "";

$_=<>;
print OUT  ("Loci\taus\tnum\tindica\tnum\taus_indica\tnum\tjaponica\tnum\tindica\tnum\tjap_indica\tnum\n");
while(<>)
{
	chomp;
	my (@e) = split;
	
	my $cmd = "";

#    my $ref = shift(@e);


	$cmd .= "$e[0]\t";
	#print OUT $chr, "\t", $loci, "\t";
	
	#for my $groups([@group_indica], [@group_japonica], [@group_nivara], [@group_rufipogon], [@group_wild_rice])

    my %base_count_raw;
    for my $i(1..$#e)
    {
        next if ($e[$i] eq 'N');
        my $base = $Abbrev{$e[$i]};#in case meeting heterozygosity
        for my $j(@$base)
        {
            $base_count_raw{$j} ++;
        }
    }
    my @best = sort {$base_count_raw{$b} > $base_count_raw{$a}} keys %base_count_raw;
    #print @best;#exit;

	for my $groups([@group_aus], [@group_indica], [@group_aus_indica], [@group_japonica], [@group_indica], [@group_jap_indica])
	{
        my $base_sequenced = 0;
        my %base_freq = (A =>0, 
        T => 0,
        C => 0,
        G => 0);
        my %base_count = %base_freq;
	
		for my $index (@$groups)
		{
			next if ($e[$index] eq "N");
   #         print($e[$index], "\n"); exit;
			my $base = $Abbrev{$e[$index]};#in case meeting heterozygosity
			
			for my $i(@$base)
			{
				$base_count{$i} ++;
                $base_sequenced++;
			}
		}

	
##		for my $index (sort keys %base_count)
##		{
##			$base_freq{$index} = $base_count{$index}/$base_sequenced;
##		}
	
		my $pi=0; #will be calculated as follows
#        for my $i ( keys %base_freq)
#        {
#            print $i, "\n";
#        }
###        my @best = sort {$base_freq{$b} > $base_freq{$a}} keys %base_freq;
#        print $best[0], "\t", $best[1], "\n";
#        print $base_freq{$best[0]}, "\t", $base_freq{$best[1]}, "\n";
#        exit;

#		for my $index (sort keys %base_freq)
#		{
#			for my $index2 (sort keys %base_freq)
#			{
#				next if ($index eq $index2);
        if($base_sequenced >= $minimum_accession && int(@best) > 1)
        {
            $base_sequenced = $base_count{$best[0]} + $base_count{$best[1]};
#            print(@best, "\n");
#            print($base_freq{$best[0]}, "\n");
#            print($base_freq{$best[1]}, "\n");
###            print $base_count{$best[0]}, "\n";
###            print $base_count{$best[1]}, "\n";
###            print $base_count{$best[1]}, "\n";
            if($base_sequenced >= $minimum_accession)
            {
                my $base_freq1 = $base_count{$best[0]} / ($base_count{$best[0]} + $base_count{$best[1]});
                my $base_freq2 = $base_count{$best[1]} / ($base_count{$best[0]} + $base_count{$best[1]});
                $base_sequenced = $base_count{$best[0]} + $base_count{$best[1]};
                $pi = $base_freq1 * $base_freq2;
            }
            #exit;
            #$pi += $base_freq{$best[0]} *$base_freq{$best[1]};
        }
#			}
#		}
		
#		if($base_sequenced < $minimum_accession)
#		{
#			$pi = $marker;
#		}
		$cmd .= "$pi\t$base_sequenced\t";
	}

	$cmd .="\n";

	print OUT $cmd;
}

