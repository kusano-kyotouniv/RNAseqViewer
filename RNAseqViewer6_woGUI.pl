# RNAseqViewer5_formac GUIぬき、キーワード検索と画像生成のみ抜粋したもの
# Tamaki, Kusano, Matsumoto, Shimada, et al., (2023投稿予定) 論文で使用したもの
# scat label 欄に volcano と書くとvolcano plot。
# 検索キーワード、volcanoかシンプルなscat plotの分け、グラフの両軸をどの列にするか等の情報はこのスクリプト内に直接書く


# 検索キーワード、volcano、使う行等の情報

my $keyword_input = 'Nicotianamine synthase';		# キーワード。図にハイライトが入る

my @scat_label;
$scat_label[0] = 'volcano';			# scat label 欄。例）１個めのグラフをvolcano plot 仕様とする場合は配列のゼロ番にvolcanoと書く。
$scat_label[1] = 'expr/cont';		# なお、volcano plot の計算とケタ合わせは別途必要。
my @scat_xfpkm;
my @scat_yfpkm;
$scat_xfpkm[0] = 2;$scat_yfpkm[0] = 3;	# それぞれの図で使う両軸をどの列にするか決める。３列目から順に０番、１番、２番、、、
$scat_xfpkm[1] = 1;$scat_yfpkm[1] = 0;
$scat_xfpkm[2] = 4;$scat_yfpkm[2] = 5;
$scat_xfpkm[3] = 1;$scat_yfpkm[3] = 4;
$scat_xfpkm[4] = 1;$scat_yfpkm[4] = 5;


# 以下プログラム

use GD;
my $ttfont = '/Library/Fonts/Arial Unicode.ttf';
my $rsem_file = 'tamakiarray_gene_RNAseqViewer2.txt';		# マイクロアレイデータのファイル
my $workidentifier = 'test';								# int(rand(100000));


# RSEM データを持つための変数
my @contig_fpkm;
my @contig_name;
my @contig_length;
my @contig_line;
my $contigs=0;
my $fpkms=0;
my @sample_name;
my $firstlinersem;

# 散布図画像を持つための変数
my $reset_scat_back=1;		# 散布図を描き直す必要が生じたとき
my $backimage;				# イメージオブジェクト。
my $saveimage_scat;
my $graphsize=900;			# キャンバスいっこ分のサイズ。拡縮するとき動的にするか？
my $scat_savefile;			# 画像をユーザー指定のファイルに保存するためのファイル名。デフォ値はundef。画像は自動的にテンポラリファイルに保存される。
		$scat_savefile ='RNAseqViewer_scatplot.png';
my $scat_showtext_limit = 5;

# 散布図を Local blast でハイライトするための変数
my @hitcontig_fpkm;		# FPKM の値はRSEMから。
my $hitcontig_resetflag = 1;
my @hitcontig_name;		# 他はLocalblastからロードする。
my @hitcontig_length;
my $hitcontig_length_max;
my $hitcontig_length_min;
my @hitcontig_score;
my $hitcontig_score_max;
my $hitcontig_score_min;
my @hitcontig_evalue;
my @hitcontig_seq;
my $hitcontigs=0;
my @hitcontig_plotx;	# プロット位置。canvas上の座標。[][] 二番目の引数がプロットの番号
my @hitcontig_ploty;

my $scats=5;
my $scat_newlabel;		# スキャッタープロットの入力用バッファ
my $scat_newx;
my $scat_newy;
my $scat_newx_buf = 'X'.$scat_newx;
my $scat_newy_buf = 'Y'.$scat_newy;
my $scat_hitextswitch = 1;
my $scat_plotcolor_switch = 0;
my $scat_zoom = 0.5;

# キーワード検索関連の変数
# my $keyword_input = 'Nicotianamine synthase';
my $keyword_error = 'OK';
my @keyword_hit_name;
my @keyword_hit_seq;
my @keyword_hit_fpkm;
my @keyword_hit_length;
my $keyword_hits = 0;
my $keyword_hit_length_max;
my $keyword_hit_length_min;
my @keyword_hit_plotx;	# プロット位置。canvas上の座標。[][] １番目の引数がグラフの番号
my @keyword_hit_ploty;

# 一枚で出力するときはこれをコメントアウト外す
#$scat_xfpkm[0] = 1;$scat_yfpkm[0] = 4;$scats=1;$scat_label[0] = '';$keyword_error='error';$scat_label[0] = 'volcano';

load_rsem();
search_keyword();			# ここをコメントアウトすればハイライトをオフにできる
draw_scat_background();
draw_scat();

sub load_rsem(){
	$fpkms = 0;
	$trinities = 0;
	$contigs = 0;
	$reset_scat_back = 1;

	open RSM, $rsem_file;
		# 最初の行が邪魔なのでスキップ。ついでにFPKM値が何サンプル分あるか調べる
		$firstlinersem = <RSM>;
		my @firstlinersem_split = split(/\t/,$firstlinersem);
		FST:for(my $i=2;$i<100;$i++){
			if($firstlinersem_split[$i] =~ /FPKM|TPM|expression/){
				$sample_name[$fpkms] = $firstlinersem_split[$i];
				$sample_name[$fpkms] =~ s/FPKM|TPM|expression//g;
				$fpkms++;
			}else{
				last FST;
			}
		}
#		$scatx_listbox -> delete(0,'end');
#		$scaty_listbox -> delete(0,'end');
#		for(my $i=0;$i<$fpkms;$i++){
#			$scatx_listbox -> insert($i,$i);
#			$scaty_listbox -> insert($i,$i);
#		}
		for(my $sc=0;$sc<$scats;$sc++){
			if($scat_xfpkm[$sc] >= $fpkms){$scat_xfpkm[$sc]=0;}
			if($scat_yfpkm[$sc] >= $fpkms){$scat_yfpkm[$sc]=0;}
		}
#		$scat_listbox -> delete(0,'end');
#		for(my $i=0;$i<$scats;$i++){$scat_listbox -> insert($i, 'x#'.$scat_xfpkm[$i].', y#'.$scat_yfpkm[$i].', '.$scat_label[$i]);}
		$reset_scat_back = 1;

		print STDERR 'RNAseq samples: '.$fpkms."\n";
		print STDERR 'Loading RSEM result...';

	while(my $line = <RSM>){	# ロード。計算の邪魔な|とかFPKM=0 とかをなんとかしておく
		chomp $line;
		my @s = split(/\t/,$line);
		$s[0] =~ s/\|/_/;							# 縦棒は計算式の邪魔になるので _ で置換する。
		$contig_line[$contigs] = $line;
		$contig_name[$contigs] = $s[0];
		$contig_length[$contigs] = $s[1];
		for(my $i=2;$i<2+$fpkms;$i++){
			if($s[$i] == 0){$s[$i] = 0.001;}		# FPKM = 0 だとlog計算できないので最低値のひとケタ下の値 0.001 を代わりに入れる
			$contig_fpkm[$contigs][$i-2] = $s[$i];
		}
		$contigs++;
		if($contigs % 10000 == 0){print STDERR '.';}
	}
	close RSM;
	print STDERR $contigs."\n";
	$rsem_file_status = $contigs. ' entries';
#	draw_scat();
}

sub search_keyword(){
	if($query_length == 0){
	if($keyword_error eq 'OK'){
		@keyword_hit_name = ();
		@keyword_hit_seq = ();			# 塩基配列を探すのはボタンを押してからにすることにした
		@keyword_hit_fpkm = ();
		@keyword_hit_length = ();
		$keyword_hits = 0;
		
		my $filename_key = $keyword_input;
		$filename_key =~ s/\|/_/g;
		$filename_key =~ s/\ /_/g;
		$filename_key =~ s/\(/_/g;
		$filename_key =~ s/\)/_/g;
		$filename_key =~ s/\,/_/g;
		open OUT, '>temp_keyhits'.$workidentifier.'_'.$filename_key.'.txt';
		print OUT $firstlinersem;
		for(my $i=0;$i<$contigs;$i++){
			if($contig_line[$i] =~ /$keyword_input/){
				print OUT $contig_line[$i]."\n";
				$keyword_hit_name[$keyword_hits] = $contig_name[$i];
				$keyword_hit_length[$keyword_hits] = $contig_length[$i];
				for(my $j=0;$j<$fpkms;$j++){$keyword_hit_fpkm[$keyword_hits][$j] = $contig_fpkm[$i][$j];}
				$keyword_hits++;
			}
		}
		close OUT;
		print STDERR 'keyword search: '.$keyword_hits." contigs hit\n";
		$keyword_error = $keyword_hits.' contigs hit';		# ヒットしたコンティぐ件数を表示する。配列の検索はユーザーが判断する

		# length が長い順に並べる
		if($keyword_hits>0){
			my @temp_name;
			my @temp_length;
			my @temp_fpkm;
			for(my $i=0;$i<$keyword_hits;$i++){
				my $max = 0;
				my $maxid;
				for(my $j=0;$j<$keyword_hits;$j++){
					if($max < $keyword_hit_length[$j]){
						$max = $keyword_hit_length[$j];
						$maxid = $j;
					}
				}
				$temp_name[$i] = $keyword_hit_name[$maxid];
				$temp_length[$i] = $keyword_hit_length[$maxid];
				for(my $k=0;$k<$fpkms;$k++){$temp_fpkm[$i][$k] = $keyword_hit_fpkm[$maxid][$k];}
				$keyword_hit_length[$maxid] = 0;
			}
			for(my $i=0;$i<$keyword_hits;$i++){
				$keyword_hit_name[$i] = $temp_name[$i];
				$keyword_hit_length[$i] = $temp_length[$i];
				for(my $k=0;$k<$fpkms;$k++){$keyword_hit_fpkm[$i][$k] = $temp_fpkm[$i][$k];}
			}
			$keyword_hit_length_max = $keyword_hit_length[0];
			$keyword_hit_length_min = $keyword_hit_length[$keyword_hits-1];
		}

		# blast結果を表示する窓に表示する
#		if($keyword_hits>0){
#			$textwidget2 -> delete('1.0','end');
#			for(my $i=0;$i<$keyword_hits;$i++){
#				$textwidget2 -> insert('end','>'.$keyword_hit_name[$i]."\n");
#				$textwidget2 -> insert('end',$keyword_hit_seq[$i]."\n");
#			}
#		}else{
#			$textwidget2 -> delete('1.0','end');
#			$textwidget2 -> insert('end','Not hits found for keyword: '.$keyword_input."\n");
#		}
#		draw_scat();
	}}
}


sub draw_scat_background(){
	if($reset_scat_back == 1){if($contigs>0){if($scats>0){

		my $z = $scat_zoom;
		$graphsize=900*$z;	# キャンバスいっこ分のサイズ
		$backimage = GD::Image ->new($graphsize * $scats,$graphsize);	# 設定数分の枠を用意する
		my $white = $backimage -> colorAllocate(255,255,255);
		my $black = $backimage -> colorAllocate(0,0,0);
		my $red = $backimage -> colorAllocate(255,0,0);
		my $green = $backimage -> colorAllocate(0,127,0);
		my $blue = $backimage -> colorAllocate(127,127,255);
		my $darkblue = $backimage -> colorAllocate(0,0,255);
		my $gray = $backimage -> colorAllocate(127,127,127);
		#$backimage -> transparent($white);
		$backimage -> interlaced('true');
		my $mergin = 80*$z;		# 左と下の余白
		my $scale = 100*$z;		# 10倍目盛のサイズ

		for(my $sc =0; $sc<$scats; $sc++){

			# 目盛線と枠を表示
			my $xzero = $mergin + $graphsize * $sc;
			my $yzero = $graphsize-$mergin;
			{
				my $fontsize = 16*$z;
				for(my $i=1;$i<8;$i++){
					$backimage -> line($xzero+$i*$scale,     $yzero,          $xzero+$i*$scale,      $yzero-8*$scale, $gray);		# 横軸補助目盛線
					$backimage -> line($xzero+$i*$scale,     $yzero,          $xzero+$i*$scale,      $yzero-8,             $black);		# 横軸めもり
					$backimage -> line($xzero+1+$i*$scale, $yzero,          $xzero+$i*$scale+1,  $yzero-8,             $black);		# 横軸めもり
					$backimage -> line($xzero,          $yzero-$i*$scale,     $xzero+8*$scale, $yzero-$i*$scale,      $gray);		# 縦軸補助目盛線
					$backimage -> line($xzero,          $yzero-$i*$scale,     $xzero+8            , $yzero-$i*$scale,      $black);		# 縦軸めもり
					$backimage -> line($xzero,          $yzero-$i*$scale-1,  $xzero+8            , $yzero-$i*$scale-1,  $black);		# 縦軸めもり
					my $label = 10**($i-3);
					if($scat_label[$sc] =~ /olcano/){$label=($i-4);}
					$backimage -> stringFT($black, $ttfont, $fontsize, 0, $xzero+4*$z+$i*$scale-length($label)*$fontsize/2,$yzero+8*$z+$fontsize,$label);		# 横軸の数字
					if($scat_label[$sc] =~ /olcano/){$label=$i;}
#					if($scat_label[$sc] =~ /olcano/){$label='1e-0'.$i;}
#					if($scat_label[$sc] =~ /olcano/){$label=1/10**$i;}
					$backimage -> stringFT($black, $ttfont, $fontsize, 0, $xzero-10*$z-length($label)*$fontsize/(1.5),$yzero-$i*$scale+$fontsize/2,$label);		# 縦軸の数字
				}
				$backimage -> line($xzero, $yzero, $xzero+8*$scale, $yzero,         $black);	# 横の目盛線
				$backimage -> line($xzero, $yzero, $xzero,          $yzero-8*$scale,$black);	# 縦の目盛線
			}
			# 軸ラベルを表示
			{	my $fontsize = 24*$z;
#				my $label = 'Sample #'.$scat_xfpkm[$sc].' (FPKM)';	# 横軸のラベル
				my $label = $sample_name[$scat_xfpkm[$sc]];			# 横軸のラベル
				   if($scat_label[$sc] =~ /olcano/){$label=$label.'(log4)';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale/2-length($label)*$fontsize/3,$yzero+$fontsize+40*$z,$label);
#				   $label = 'Sample #'.$scat_yfpkm[$sc].' (FPKM)';	# 縦軸のラベル
				   $label = $sample_name[$scat_yfpkm[$sc]];			# 縦軸のラベル
				   if($scat_label[$sc] =~ /olcano/){$label=$label.'(-log10)';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 1.5708, $xzero-$fontsize-32*$z, $yzero-8*$scale/2+length($label)*$fontsize/3,$label);
				   $label = $scat_label[$sc];						# 図のラベル
				   if($scat_label[$sc] =~ /olcano/){$label='';}
				$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale-length($label)*$fontsize/1.5,$yzero-$fontsize,$label);
			}

			# 背景をプロットする
			my $pointsize = 2;
			for(my $i=0;$i<$contigs;$i++){
				my $x1 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
				my $y1 = $yzero - $scale * (log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
				my $x2 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
				my $y2 = $yzero - $scale * (log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
				if($scat_label[$sc] =~ /olcano/){
					$x1 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
					$y1 = $yzero - $scale * (-log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
					$x2 = $xzero + $scale * (log($contig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
					$y2 = $yzero - $scale * (-log($contig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
				}
				$backimage -> filledRectangle($x1,$y1,$x2,$y2,$blue);
			}
			
			# x=y の線を引く
			if($scat_label[$sc] =~ /olcano/){
				$backimage -> line($xzero+4*$scale,   $yzero,   $xzero+4*$scale,   $yzero-8*$scale,   $green);
				$backimage -> line($xzero+4*$scale+1, $yzero,   $xzero+4*$scale+1, $yzero-8*$scale,   $green);	# ちょっと細いので太くする
				$backimage -> line($xzero+4*$scale-1, $yzero,   $xzero+4*$scale-1, $yzero-8*$scale,   $green);
			}else{
				$backimage -> line($xzero,   $yzero,   $xzero+8*$scale,   $yzero-8*$scale,   $green);
				$backimage -> line($xzero+1, $yzero,   $xzero+8*$scale,   $yzero-8*$scale+1, $green);	# ちょっと細いので太くする
				$backimage -> line($xzero,   $yzero-1, $xzero+8*$scale-1, $yzero-8*$scale,   $green);
			}
		}
		$reset_scat_back = 0;
	}}}
}


sub draw_scat(){
	print STDERR 'call draw_scat'."\n";

			my $z = $scat_zoom;
			$graphsize=900*$z;	# キャンバスいっこ分のサイズ
			my $mergin = 80*$z;	# 左と下の余白
			my $scale = 100*$z;	# 10倍目盛のサイズ
			my $pointsize = 5;
			my @fieldx;
			my @fieldy;
			my $fields=0;

	if($contigs>0){if($scats>0){
		if($reset_scat_back == 1){draw_scat_background();}			# 背景画像のリセットフラグの参照が二重だが一応とっとく
#		if($hitcontig_resetflag == 1){search_hitcontig_fpkm();}		# Local blast が更新された場合は描画前にfpkm値の検索をやり直す
#		if($textwidget_refreshflag == 1){input_fasta();}			# スキャッタープロットの前にエントリーリストを確認する。
#		if(Exists($plot_window)){}else{
#			$plot_window = $control_window -> Toplevel();
#			$plot_window -> title('ScatterPlot View '.$workidentifier);
#			$plot_window -> geometry("1300x650");
#			$plot_frame = $plot_window -> Scrolled('Canvas',-scrollbars => 'se') -> pack(-fill => 'both', -expand => 1);
#			$plot_canvas = $plot_frame -> Subwidget('canvas') -> pack(-fill => 'both', -expand => 1);
#		}
#		my $png = $plot_window -> Photo(-data => encode_base64($backimage -> png), -format => 'png'); # Photo が機能しないので
#		$plot_canvas -> delete('all');
#		$plot_frame -> configure (-width => $graphsize * $scats, -height => $graphsize);
#		$plot_canvas -> configure (-width => $graphsize * $scats, -height => $graphsize);	# これだけでは何故かサイズを小さい方には変更できない。
#		$plot_canvas -> configure (-scrollregion => [0,0,$graphsize * $scats,$graphsize]);
#		my $background_plotimage = $plot_canvas -> create('image', $graphsize * $scats/2, $graphsize/2, -image => $png);
		# GDの像を表示できないので代わりに枠線と斜め線を引く
#		for(my $sc=0;$sc<$scats;$sc++){
#			my $xzero = $mergin + $graphsize * $sc;
#			my $yzero = $graphsize-$mergin;
#			$plot_canvas -> create('rectangle', $xzero, $yzero, $xzero+8*$scale, $yzero-8*$scale, -fill => 'white', -outline => 'black',-width => 1);
#				my $fontsize = 16*$z;
#				for(my $i=1;$i<8;$i++){
#					$plot_canvas -> create('line', $xzero+$i*$scale,$yzero,          $xzero+$i*$scale,$yzero-8*$scale, -fill => 'gray', -width => 1);		# 横軸補助目盛線
#					$plot_canvas -> create('line', $xzero+$i*$scale,$yzero,          $xzero+$i*$scale,$yzero-8,        -fill => 'black', -width => 2);		# 横軸めもり
#					$plot_canvas -> create('line', $xzero,          $yzero-$i*$scale,$xzero+8*$scale, $yzero-$i*$scale, -fill => 'gray', -width => 1);		# 縦軸補助目盛線
#					$plot_canvas -> create('line', $xzero,          $yzero-$i*$scale,$xzero+8       , $yzero-$i*$scale, -fill => 'black', -width => 2);		# 縦軸めもり
#					my $label = 10**($i-3);
#					$plot_canvas -> create('text', $xzero+$i*$scale-length($label)*$fontsize/16,$yzero+$fontsize, -text => $label, -font => 'Arial '.$fontsize); # 横軸の数字
#					$plot_canvas -> create('text', $xzero-14*$z-length($label)*$fontsize/4,$yzero-$i*$scale, -text => $label, -font => 'Arial '.$fontsize); # 縦軸の数字
#				}
#			$plot_canvas -> create('line', $xzero, $yzero, $xzero+8*$scale, $yzero-8*$scale, -fill => 'blue', -width => 2);
#			$fontsize = 24*$z;
#			my $label = 'Sample #'.$scat_yfpkm[$sc].' / #'.$scat_xfpkm[$sc].': '.$scat_label[$sc];						# 図のラベル
#			$plot_canvas -> create('text', $xzero+7*$scale-length($label)*$fontsize/4,$yzero-$fontsize-$scale*0.3, -text => $label, -font => 'Arial '.$fontsize);
#			$backimage -> stringFT($black, $ttfont, $fontsize, 0,      $xzero+8*$scale-length($label)*$fontsize/1.5,$yzero-$fontsize,$label);
#		}
#		if($query_length == 0 && $keyword_hits == 0){$plot_canvas -> bind($background_plotimage, "<Button-1>" => [\&create_by_userclick, Ev('x'), Ev('y')]);}
#		$plot_window -> maxsize($graphsize * $scats+20,$graphsize+20);	# これで一応小さい方にも変更できるが、一回ウィンドウに触ってもらう必要がある。仕方ないからこれでいい

		$saveimage_scat = $backimage -> clone();
		my $white = $saveimage_scat -> colorAllocate(255,255,255);
		my $black = $saveimage_scat -> colorAllocate(0,0,0);
		my $red = $saveimage_scat -> colorAllocate(255,0,0);
		my $savecolor_index;
		my @savecolor;
		for(my $i=0;$i<128;$i++){
			my $r =     int(511* $i /127);if($r>255){$r=255;}
			my $g = 511-int(511* $i /127);if($g>255){$g=255;}
			my $b = 0;
			$savecolor[$i] = $saveimage_scat -> colorAllocate($r,$g,$b);
		}


		# Local Blast の結果をプロットに反映させる

		if($query_length > 1){
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# プロットする
				for(my $i=$hitcontigs-1;$i>=0;$i--){
					my $x1 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
					my $y1 = $yzero - $scale * (log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
					my $x2 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
					my $y2 = $yzero - $scale * (log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
					if($scat_label[$sc] =~ /olcano/){
						$x1 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
						$y1 = $yzero - $scale * (-log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
						$x2 = $xzero + $scale * (log($hitcontig_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
						$y2 = $yzero - $scale * (-log($hitcontig_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
					}

					my $pointcolor;
					if($scat_plotcolor_switch == 0){
						my $r = int(511 * ($hitcontig_score[$i] - $hitcontig_score_min)/ ($hitcontig_score_max - $hitcontig_score_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($hitcontig_score[$i] - $hitcontig_score_min) / ($hitcontig_score_max - $hitcontig_score_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127*($hitcontig_score[$i] - $hitcontig_score_min)/ ($hitcontig_score_max - $hitcontig_score_min));
					}
					if($scat_plotcolor_switch == 1){
						my $r = int(511 * $hitcontig_evalue[$i] / 181);
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * $hitcontig_evalue[$i] / 181);
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * $hitcontig_evalue[$i] / 181);
					}
					if($scat_plotcolor_switch == 2){
						my $r = int(511 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($hitcontig_length[$i] - $hitcontig_length_min) / ($hitcontig_length_max - $hitcontig_length_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
					}
#					my $plot = $plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -fill => $pointcolor, -width => 2);
					$saveimage_scat -> filledRectangle($x1,$y1,$x2,$y2,$black);
					$saveimage_scat -> filledRectangle($x1+1,$y1+1,$x2-1,$y2-1,$savecolor[$savecolor_index]);
#					$plot_canvas -> bind($plot, "<Button-1>" => [\&select_by_userclick, Ev('x'), Ev('y')]);
					$hitcontig_plotx[$sc][$i] = ($x1+$x2)/2;	# find が機能しないので仕方なく
					$hitcontig_ploty[$sc][$i] = ($y1+$y2)/2;

					# プロットした座標をとっとく
					$fieldx[$fields]=($x1+$x2)/2;
					$fieldy[$fields]=($y1+$y2)/2;
					$fields++;
				}
			}
		}

		if($query_length == 0){if($keyword_hits>0){		# キーワード検索のときの色
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# プロットする
				for(my $i=$keyword_hits-1;$i>=0;$i--){
					my $x1 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
					my $y1 = $yzero - $scale * (log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
					my $x2 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
					my $y2 = $yzero - $scale * (log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
					if($scat_label[$sc] =~ /olcano/){
						$x1 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
						$y1 = $yzero - $scale * (-log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
						$x2 = $xzero + $scale * (log($keyword_hit_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
						$y2 = $yzero - $scale * (-log($keyword_hit_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
					}

					my $pointcolor;
					if($keyword_hit_length_max - $keyword_hit_length_min > 0){		# ゼロ割エラー出ちゃうので。
						my $r = int(511 * ($keyword_hit_length[$i] - $keyword_hit_length_min)/ ($keyword_hit_length_max - $keyword_hit_length_min));
						if($r>255){$r = 255;};$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g = 511-int(511 * ($keyword_hit_length[$i] - $keyword_hit_length_min) / ($keyword_hit_length_max - $keyword_hit_length_min));
						if($g>255){$g = 255;};$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b = sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = int(127 * ($keyword_hit_length[$i] - $keyword_hit_length_min)/ ($keyword_hit_length_max - $keyword_hit_length_min));
#						$savecolor_index = int(127 * ($hitcontig_length[$i] - $hitcontig_length_min)/ ($hitcontig_length_max - $hitcontig_length_min));
					}else{
						my $r=255;$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g=255;$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b=sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = 63;
					}

					if($scat_plotcolor_switch == 0){	# キーワード検索のハイライトの場合に黄色一色にするスイッチとして使用してみる
						my $r=255;$r=sprintf("%X",$r);if(length($r)==1){$r='0'.$r;}
						my $g=255;$g=sprintf("%X",$g);if(length($g)==1){$g='0'.$g;}
						my $b = 0;$b=sprintf("%X",$b);if(length($b)==1){$b='0'.$b;}
						$pointcolor = '#'.$r.$g.$b;
						$savecolor_index = 63;
					}

#					my $plot = $plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -fill => $pointcolor, -width => 2);
					$saveimage_scat -> filledRectangle($x1,$y1,$x2,$y2,$black);
					$saveimage_scat -> filledRectangle($x1+2,$y1+2,$x2-2,$y2-2,$savecolor[$savecolor_index]);
#					$plot_canvas -> bind($plot, "<Button-1>" => [\&select_by_userclick, Ev('x'), Ev('y')]);
					$keyword_hit_plotx[$sc][$i] = ($x1+$x2)/2;	# find が機能しないので仕方なく
					$keyword_hit_ploty[$sc][$i] = ($y1+$y2)/2;

					# プロットした座標をとっとく
					$fieldx[$fields]=($x1+$x2)/2;
					$fieldy[$fields]=($y1+$y2)/2;
					$fields++;
				}
				
				# キーワードを表示する
				my $label = $keyword_input.' ('.$keyword_hits.' genes)';
				$label =~ s/\|/\,\ /g;			# or 演算子は検索に有効だけど見づらいのでコンマと空白を入れる
				my $title_x = $xzero+$scale*1;	# 左から１マス目
				my $title_y = $yzero-$scale*6;	# 下から６マス目
				my $fontsize = 24*$z;			# 軸ラベルと同じサイズ= 24 * $z
				$saveimage_scat -> stringFT($black, $ttfont, $fontsize, 0, $title_x, $title_y, $label);
			}
		}}

		if($fastas > 0){
			for(my $sc =0; $sc<$scats; $sc++){
				my $xzero = $mergin + $graphsize * $sc;
				my $yzero = $graphsize-$mergin;
				# ユーザーの指定コンティグをハイライト表示にする
					for(my $i=$fastas-1;$i>=0;$i--){
						if(defined($fasta_fpkm[$i][$scat_xfpkm[$sc]])){
							my $x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
							my $y1 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
							my $x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
							my $y2 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
							if($scat_label[$sc] =~ /olcano/){
								$x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
								$y1 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
								$x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
								$y2 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
							}
#							$plot_canvas -> create('rectangle', $x1,  $y1,  $x2,  $y2,  -outline => 'black', -width => 2);
							$saveimage_scat -> rectangle($x1,$y1,$x2,$y2,$black);
							$saveimage_scat -> rectangle($x1+2,$y1+2,$x2-2,$y2-2,$black);
							$x1 = $x1 - $pointsize;$y1 = $y1 - $pointsize;
							$x2 = $x2 + $pointsize;$y2 = $y2 + $pointsize;
#							$plot_canvas -> create('oval',      $x1,  $y1,  $x2,  $y2,  -outline => 'red', -width => 2);
							$saveimage_scat -> arc(($x1+$x2)/2, ($y1+$y2)/2, $pointsize*4, $pointsize*4, 0, 360, $red);
							$fieldx[$fields]=($x1+$x2)/2;
							$fieldy[$fields]=($y1+$y2)/2;
							$fields++;
						}
					}
				my $centerx=0;
				my $centery=0;
				my $centerc=0;
				if($fields>0){
					# 文字をずらす方向を決めるために今あるスポット中心あたりの座標を探す
					for(my $i=0;$i<$fields;$i++){
						if($fieldx[$i]>$xzero){if($fieldy[$i]<$yzero){
							$centerx = $centerx + $fieldx[$i];
							$centery = $centery + $fieldy[$i];
							$centerc++;
						}}
					}
					if($centerc>0){
						$centerx = $centerx / $centerc;
						$centery = $centery / $centerc;
					}else{
						$centerx = $xzero + $scale * 4;
						$centery = $yzero - $scale * 4;
					}
				}else{
					$centerx = $xzero + $scale * 4;
					$centery = $yzero - $scale * 4;
				}
					for(my $i=$fastas-1;$i>=0;$i--){
						if(defined($fasta_fpkm[$i][$scat_xfpkm[$sc]])){
							if($scat_hitextswitch == 1){
							if($i<$scat_showtext_limit){

								# 文字が目立つ位置まで移動する
								my $x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)-$pointsize;
								my $y1 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)-$pointsize;
								my $x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(10)+3)+$pointsize;
								my $y2 = $yzero - $scale * (log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+3)+$pointsize;
								if($scat_label[$sc] =~ /olcano/){
									$x1 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)-$pointsize;
									$y1 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)-$pointsize;
									$x2 = $xzero + $scale * (log($fasta_fpkm[$i][$scat_xfpkm[$sc]])/log(4)+4)+$pointsize;
									$y2 = $yzero - $scale * (-log($fasta_fpkm[$i][$scat_yfpkm[$sc]])/log(10)+0)+$pointsize;
								}
								my $posx = ($x1+$x2)/2;
								my $posy = ($y1+$y2)/2;
								my $effort =0;
								my $nearest=0;
								DST:while($nearest<150*$z){
									# 良い感じのコンティグの場合
									if($fasta_fpkm[$i][$scat_xfpkm[$sc]] >0.001 && $fasta_fpkm[$i][$scat_yfpkm[$sc]] >0.001){
										my $vector = sqrt(($posx - $centerx)**2 + ($posy - $centery)**2);	# スポット群の中心点から遠ざかる方向へのベクトルの長さ
										my $force = 0;
										if($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]] != 0){			# 右下か左上へのベクトル
											$force = ($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]])
											    / abs($fasta_fpkm[$i][$scat_xfpkm[$sc]] - $fasta_fpkm[$i][$scat_yfpkm[$sc]]);
										}
										$posx = $posx + (10 * ($posx - $centerx)/$vector + rand(40)-20 + 6*$force)*$z;		# x方向のランダムを気持ち多め
										$posy = $posy + (10 * ($posy - $centery)/$vector + rand(30)-15 + 6*$force)*$z;
										if($posx < $xzero + 80*$z){$posx = ($x1+$x2)/2;$posy = ($y1+$y2)/2;}				# ハミ出たら戻す
										if($posy > $yzero - 12*$z){$posy = ($y1+$y2)/2;$posx = ($x1+$x2)/2;}
										if($posx > $xzero + $scale * 8 - 80*$z){$posx = $xzero + $scale * 8 - 80*$z;}		# ハミ出たら戻す
										if($posy < $yzero - $scale * 8 + 12*$z){$posy = $yzero - $scale * 8 + 12*$z;}
									# どっちかの軸がゼロになっている変なコンティグの場合
									}else{
										if($fasta_fpkm[$i][$scat_xfpkm[$sc]] == 0.001 && $fasta_fpkm[$i][$scat_yfpkm[$sc]] == 0.001){
											$posx = $posx + rand(10)*$z;$posy = $posy - rand(10)*$z;	# 両軸ともゼロの場合は右上方向に移動
										}else{
											if($fasta_fpkm[$i][$scat_xfpkm[$sc]] == 0.001){
												$posx = $posx + rand(10)*$z;						# x が初めからゼロの場合は右方向に移動
												$posy = $posy + (rand(20)-10)*$z;
												if($posx > $xzero + $scale * 2 - 80*$z){$posx = ($x1+$x2)/2;}		# ハミ出たら戻す
											}
											if($fasta_fpkm[$i][$scat_yfpkm[$sc]] == 0.001){
												$posx = $posx + (rand(20)-10)*$z;
												$posy = $posy - rand(10)*$z;
												if($posy < $yzero - $scale * 1 + 12*$z){$posy = ($y1+$y2)/2;}
											}
										}

									}

									# 一番近いプロットまでの距離を雑に求めて評価用の変数に入れる
									my $dist_min=1000;
									for(my $j=0;$j<$fields;$j++){
										my $dist = 80*$z+sqrt(abs($posx - $fieldx[$j])**2 + abs($posy - $fieldy[$j])**2);
										if($dist_min > $dist){$dist_min = $dist;}
									}
									$nearest = $dist_min;
									$effort++;
									if($effort>1000){last DST;}
								}
								my $text = $fasta_name[$i]; $text =~ s/_//g;
#								$plot_canvas -> create('line', $posx,$posy, ($x1+$x2)/2, ($y1+$y2)/2, -width => 2, -fill => 'red');
					#			$plot_canvas -> create('rectangle', $posx-80*$z,$posy-12*$z,$posx+80*$z,$posy+12*$z, -width => 2, -fill => 'white', -outline => 'red');
#								$plot_canvas -> create('rectangle', $posx-80,$posy-12,$posx+80,$posy+12, -width => 2, -fill => 'white', -outline => 'red');
								$saveimage_scat -> line($posx,$posy, ($x1+$x2)/2, ($y1+$y2)/2,$red);
					#			$saveimage_scat -> filledRectangle($posx-80*$z-1,$posy-12*$z-1,$posx+80*$z+1,$posy+12*$z+1,$red);
					#			$saveimage_scat -> filledRectangle($posx-80*$z+1,$posy-12*$z+1,$posx+80*$z-1,$posy+12*$z-1,$white);
								$saveimage_scat -> filledRectangle($posx-80-1,$posy-12-1,$posx+80+1,$posy+12+1,$red);
								$saveimage_scat -> filledRectangle($posx-80+1,$posy-12+1,$posx+80-1,$posy+12-1,$white);
					#			my $fontsize = 14*$z;
								my $fontsize = 14;
#								$plot_canvas -> create('text',$posx,$posy, -text => $text, -font => 'Arial '.$fontsize);
								$saveimage_scat -> stringFT($black, $ttfont, $fontsize, 0, $posx-length($text)*$fontsize/3,$posy+$fontsize/2, $text);
								$fieldx[$fields]=$posx;
								$fieldy[$fields]=$posy;
								$fields++;
							}}
						}
					}
				#}
			}
		}

		# セーブボタンを表示
#		my $plot_savebutton_rec = $plot_canvas -> create('rectangle', 4,4,4+48,4+24, -outline => 'black', -fill => 'gray');
#		   $plot_canvas -> bind($plot_savebutton_rec, "<Button-1>" => [\&save_scat]);
#		my $plot_savebutton_txt = $plot_canvas -> create('text', (4+48)/2,(4+24)/2, -anchor => 'c', -text => 'Save');
#		   $plot_canvas -> bind($plot_savebutton_txt, "<Button-1>" => [\&save_scat]);

#		open OUT, ">temp_drawplot.png";
#		binmode OUT;
#		print OUT $saveimage_scat -> png();
#		close OUT;
#
		if(defined($scat_savefile)){
			if($scat_savefile =~ /\.png/){}else{$scat_savefile = $scat_savefile.'.png';}
			open OUT, '>'.$scat_savefile;
			binmode OUT;
			print OUT $saveimage_scat -> png();
			close OUT;
			print STDERR 'save done: '.$scat_savefile."\n";
#			undef($scat_savefile);
		}

	}}
}

