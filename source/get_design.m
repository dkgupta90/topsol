function [xleaf, nely, nelx, fixeddofs, nodofs] = get_design(design)
    if design == 0
        xleaf = zeros(401);
        xleaf(199:201, 199:201) = 100;
        nely = 401 - 1;
        nelx = 401 - 1;
    elseif design == 1
        xleaf = load('leaf/leaf1.dat', '-ascii');
        xleaf = xleaf(1:860, 225:650);
        xleaf(859:860, 210:225) = 100;
        %xleaf(835:860, 195:230) = 100;
        xleaf(750:800, 210:240) = 0;
        nely = 860 - 1;
        nelx = 426 - 1;
    elseif design == 2
        xleaf = load('leaf/leaf_maple2.dat', '-ascii');
        xleaf(348:350, 302:304) = 100;
        nely = 726 - 1;
        nelx = 628 - 1; %this is wrong, correct it
    elseif design == 3   %%%Hexagon
        xleaf = load('leaf/hexagon1.dat', '-ascii');
        xleaf(314:316, 2:4) = 100;   %node1
        nely = 628 - 1;
        nelx = 726 - 1;
    elseif design == 4   %%%Hexagon
        xleaf = load('leaf/hexagon1.dat', '-ascii');
        xleaf(314:316, 2:4) = 100;   %node1
        xleaf(1:3, 182:184) = 100;   %node2
        xleaf(1:3, 543:545) = 100;   %node3
        xleaf(626:628, 182:184) = 100;   %node4
        xleaf(626:628, 543:545) = 100;   %node5    
        xleaf(314:316, 723:725) = 100; %node6
        nely = 628 - 1;
        nelx = 726 - 1;
    elseif design == 5
        xleaf = load('leaf/bike_head.dat', '-ascii');
        xleaf(xleaf < 0) = 0;
        %xleaf = 255 - xleaf;
        xleaf(40:51, 505:512) = 0;
        xleaf(40:51, 512:525) = 255;
        xleaf(360:390, 550:570) = 0;
        %corrections above
        xleaf(460:462, 320:329) = 100;
        nely = 495 - 1;
        nelx = 636 - 1;
    elseif design == 6
        xleaf = load('leaf/bike_head.dat', '-ascii');
        xleaf(xleaf < 0) = 0;
        %xleaf = 255 - xleaf;
        xleaf(40:51, 505:512) = 0;
        xleaf(40:51, 512:525) = 255;
        xleaf(360:390, 550:570) = 0;
        %corrections above
        xleaf(460:462, 320:329) = 100;
        xleaf(25:27, 130:132) = 100;
        xleaf(21:23, 490:492) = 100;
        xleaf(130:132, 250:252) = 100;
        xleaf(130:132, 375:377) = 100;
        xleaf(200:202, 100:102) = 100;
        xleaf(200:202, 533:535) = 100;
        xleaf(475:477, 110:112) = 100;
        xleaf(463:465, 525:527) = 100;
        xleaf(110:112, 150:152) = 100;
        xleaf(110:112, 480:482) = 100;
        xleaf(400:402, 230:232) = 100;
        xleaf(400:402, 410:412) = 100;
        xleaf(186:188, 318:322) = 100;
        xleaf(322:324, 52:54) = 100;
        xleaf(308:310, 580:582) = 100;
        nely = 495 - 1;
        nelx = 636 - 1;
    elseif design == 7
        xleaf = load('leaf/bike_side.dat', '-ascii');
        xleaf = kron(xleaf, ones(2, 2));
        xleaf(xleaf < 0) = 0;
        xleaf = 255 - xleaf;
        xleaf(390:392, 24:26) = 100;
        nely = 422 - 1;
        nelx = 672 - 1;
    elseif design == 8
        xleaf = load('leaf/bike_side.dat', '-ascii');
        xleaf = kron(xleaf, ones(2, 2));
        xleaf(xleaf < 0) = 0;
        xleaf = 255 - xleaf;
        xleaf(390:392, 24:26) = 100;
        xleaf(80:82, 660:662) = 100;
        xleaf(395:397, 505:507) = 100;
        xleaf(350:352, 150:152) = 100;
        xleaf(300:302, 300:302) = 100;
        xleaf(200:202, 400:402) = 100;
        xleaf(120:122, 350:352) = 100;
        xleaf(80:82, 497:499) = 100;
        xleaf(125:127, 225:227) = 100;
        xleaf(195:197, 288:290) = 100;
        xleaf(383:385, 248:250) = 100;
        xleaf(325:327, 379:381) = 100;
        xleaf(160:162, 560:562) = 100;
        nely = 422 - 1;
        nelx = 672 - 1;
    elseif design == 9
        xleaf = load('leaf/circle.dat', '-ascii');
        xleaf = kron(xleaf, ones(2, 2));
        xleaf(207:210, 207:209) = 100;
        nely = 416 - 1;
        nelx = 414 - 1;
    elseif design == 10
        xleaf = imread('leaf/christmas_tree/Christmas_Tree_v3-adjustedcolours2.png');
        xleaf = im2double(xleaf);
        xleaf(xleaf == 1) = 255;
        xleaf(:, 505:-1:254) = xleaf(:, 1:252);
        xleaf(10:11, 249:255) = 100;
        nely = 661 - 1;
        nelx = 505 - 1;
%         xleaf(1320:1321, 494:515) = 100;
%         nely = 1321 - 1;
%         nelx = 1009 - 1;
    elseif design == 11
        % This is for the c-Si hexagonal shape solar cell
        nelx = 400;
        nely = 400;
        [freedofs, fixeddofs] = generate_hexcell_freedofs(nelx+1, nely+1, 40);
        nodofs = setdiff(1:(nelx+1)*(nely+1), freedofs);
        nodofs = setdiff(nodofs, fixeddofs);

        xleaf = ones(1, (nelx+1)*(nely+1));
        xleaf(fixeddofs) = 100;
        xleaf(nodofs) = 255;
        xleaf = reshape(xleaf, nely+1, nelx+1);
    end
    imagesc(xleaf);
    fixeddofs = find(xleaf == 100);
    nodofs = find(xleaf == 255);
end