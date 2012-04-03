function exc = usd_btc_virwox(sll_price, btc_price, d_s_fee, s_b_fee, to_or_from)
% Gives you the Virwox exchange rate from USD to BTC including transaction
% fees. d_s_fee is the fee for getting SLL from dollars, s_b_fee is the fee
% for getting BTC from SLL.
%
% If to_or_from = 1, then you get the price of a BTC in dollars
% If to_or_from = 2, then you get the price of a USD in bitcoins

if nargin < 5
    to_or_from = 1;
end

if to_or_from == 1
    num_sll = sll_price*(1 - d_s_fee/100);
    exc = (num_sll - num_sll*s_b_fee/100)/btc_price;
end

if to_or_from == 2
    num_sll = (1 - s_b_fee/100)*btc_price;
    exc = (num_sll - num_sll*d_s_fee/100)/sll_price;
end