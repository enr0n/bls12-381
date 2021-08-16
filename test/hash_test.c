#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "finite_field.h"

#include "hash.h"

#include "BLS12_381.h"

#define TEST_MAIN_INIT                  \
    int result, pass_count, fail_count; \
                                        \
    pass_count = fail_count = 0;

#define TEST_MAIN_RETURN                                                         \
    printf("Hash test results: %d passed, %d failed\n", pass_count, fail_count); \
                                                                                 \
    return !!(fail_count);

#define TEST_RUN(test)                   \
    result = test;                       \
    if (!result) {                       \
        printf("FAIL\n\n");              \
        fail_count++;                    \
    } else {                             \
        printf("PASS\n\n");              \
        pass_count++;                    \
    }

/* Test vectors from https://datatracker.ietf.org/doc/html/draft-irtf-cfrg-hash-to-curve-11#appendix-K.1 */
#define EXPAND_MESSAGE_XMD_TEST_DST "QUUX-V01-CS02-with-expander"
#define G1_TEST_DST "QUUX-V01-CS02-with-BLS12381G1_XMD:SHA-256_SSWU_RO_"
#define G2_TEST_DST "QUUX-V01-CS02-with-BLS12381G2_XMD:SHA-256_SSWU_RO_"

const char *test_vec_msg[10] = {
    "",
    "abc",
    "abcdef0123456789",
    "q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq",
    "a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
    "",
    "abc",
    "abcdef0123456789",
    "q128_qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq",
    "a512_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
};

const uint32_t test_vec_len_in_bytes[10] = {
    0x20,
    0x20,
    0x20,
    0x20,
    0x20,
    0x80,
    0x80,
    0x80,
    0x80,
    0x80
};

const char *test_vec_uniform_bytes[10] = {
    "f659819a6473c1835b25ea59e3d38914c98b374f0970b7e4c92181df928fca88",
    "1c38f7c211ef233367b2420d04798fa4698080a8901021a795a1151775fe4da7",
    "8f7e7b66791f0da0dbb5ec7c22ec637f79758c0a48170bfb7c4611bd304ece89",
    "72d5aa5ec810370d1f0013c0df2f1d65699494ee2a39f72e1716b1b964e1c642",
    "3b8e704fc48336aca4c2a12195b720882f2162a4b7b13a9c350db46f429b771b",
    "8bcffd1a3cae24cf9cd7ab85628fd111bb17e3739d3b53f89580d217aa79526f1708354a76a402d3569d6a9d19ef3de4d0b991e4f54b9f20dcde9b95a66824cbdf6c1a963a1913d43fd7ac443a02fc5d9d8d77e2071b86ab114a9f34150954a7531da568a1ea8c760861c0cde2005afc2c114042ee7b5848f5303f0611cf297f",
    "fe994ec51bdaa821598047b3121c149b364b178606d5e72bfbb713933acc29c186f316baecf7ea22212f2496ef3f785a27e84a40d8b299cec56032763eceeff4c61bd1fe65ed81decafff4a31d0198619c0aa0c6c51fca15520789925e813dcfd318b542f8799441271f4db9ee3b8092a7a2e8d5b75b73e28fb1ab6b4573c192",
    "c9ec7941811b1e19ce98e21db28d22259354d4d0643e301175e2f474e030d32694e9dd5520dde93f3600d8edad94e5c364903088a7228cc9eff685d7eaac50d5a5a8229d083b51de4ccc3733917f4b9535a819b445814890b7029b5de805bf62b33a4dc7e24acdf2c924e9fe50d55a6b832c8c84c7f82474b34e48c6d43867be",
    "48e256ddba722053ba462b2b93351fc966026e6d6db493189798181c5f3feea377b5a6f1d8368d7453faef715f9aecb078cd402cbd548c0e179c4ed1e4c7e5b048e0a39d31817b5b24f50db58bb3720fe96ba53db947842120a068816ac05c159bb5266c63658b4f000cbf87b1209a225def8ef1dca917bcda79a1e42acd8069",
    "396962db47f749ec3b5042ce2452b619607f27fd3939ece2746a7614fb83a1d097f554df3927b084e55de92c7871430d6b95c2a13896d8a33bc48587b1f66d21b128a1a8240d5b0c26dfe795a1a842a0807bb148b77c2ef82ed4b6c9f7fcb732e7f94466c8b51e52bf378fba044a31f5cb44583a892f5969dcd73b3fa128816e"
};

const char *test_vec_G1_field_elems_0[5] = {
    "0x0ba14bd907ad64a016293ee7c2d276b8eae71f25a4b941eece7b0d89f17f75cb3ae5438a614fb61d6835ad59f29c564f",
    "0x0d921c33f2bad966478a03ca35d05719bdf92d347557ea166e5bba579eea9b83e9afa5c088573c2281410369fbd32951",
    "0x062d1865eb80ebfa73dcfc45db1ad4266b9f3a93219976a3790ab8d52d3e5f1e62f3b01795e36834b17b70e7b76246d4",
    "0x010476f6a060453c0b1ad0b628f3e57c23039ee16eea5e71bb87c3b5419b1255dc0e5883322e563b84a29543823c0e86",
    "0x0a8ffa7447f6be1c5a2ea4b959c9454b431e29ccc0802bc052413a9c5b4f9aac67a93431bd480d15be1e057c8a08e8c6"
};

const char *test_vec_G1_field_elems_1[5] = {
    "0x019b9bd7979f12657976de2884c7cce192b82c177c80e0ec604436a7f538d231552f0d96d9f7babe5fa3b19b3ff25ac9",
    "0x003574a00b109ada2f26a37a91f9d1e740dffd8d69ec0c35e1e9f4652c7dba61123e9dd2e76c655d956e2b3462611139",
    "0x0cdc3e2f271f29c4ff75020857ce6c5d36008c9b48385ea2f2bf6f96f428a3deb798aa033cd482d1cdc8b30178b08e3a",
    "0x0b1a912064fb0554b180e07af7e787f1f883a0470759c03c1b6509eb8ce980d1670305ae7b928226bb58fdc0a419f46e",
    "0x05d487032f602c90fa7625dbafe0f4a49ef4a6b0b33d7bb349ff4cf5410d297fd6241876e3e77b651cfc8191e40a68b7"
};

const char *test_vec_G2_field_elems_0_0[5] = {
    "0x03dbc2cce174e91ba93cbb08f26b917f98194a2ea08d1cce75b2b9cc9f21689d80bd79b594a613d0a68eb807dfdc1cf8",
    "0x15f7c0aa8f6b296ab5ff9c2c7581ade64f4ee6f1bf18f55179ff44a2cf355fa53dd2a2158c5ecb17d7c52f63e7195771",
    "0x0313d9325081b415bfd4e5364efaef392ecf69b087496973b229303e1816d2080971470f7da112c4eb43053130b785e1",
    "0x025820cefc7d06fd38de7d8e370e0da8a52498be9b53cba9927b2ef5c6de1e12e12f188bbc7bc923864883c57e49e253",
    "0x190b513da3e66fc9a3587b78c76d1d132b1152174d0b83e3c1114066392579a45824c5fa17649ab89299ddd4bda54935"
};

const char *test_vec_G2_field_elems_0_1[5] = {
    "0x05a2acec64114845711a54199ea339abd125ba38253b70a92c876df10598bd1986b739cad67961eb94f7076511b3b39a",
    "0x01c8067bf4c0ba709aa8b9abc3d1cef589a4758e09ef53732d670fd8739a7274e111ba2fcaa71b3d33df2a3a0c8529dd",
    "0x062f84cb21ed89406890c051a0e8b9cf6c575cf6e8e18ecf63ba86826b0ae02548d83b483b79e48512b82a6c0686df8f",
    "0x034147b77ce337a52e5948f66db0bab47a8d038e712123bb381899b6ab5ad20f02805601e6104c29df18c254b8618c7b",
    "0x12ab625b0fe0ebd1367fe9fac57bb1168891846039b4216b9d94007b674de2d79126870e88aeef54b2ec717a887dcf39"
};

const char *test_vec_G2_field_elems_1_0[5] = {
    "0x02f99798e8a5acdeed60d7e18e9120521ba1f47ec090984662846bc825de191b5b7641148c0dbc237726a334473eee94",
    "0x187111d5e088b6b9acfdfad078c4dacf72dcd17ca17c82be35e79f8c372a693f60a033b461d81b025864a0ad051a06e4",
    "0x1739123845406baa7be5c5dc74492051b6d42504de008c635f3535bb831d478a341420e67dcc7b46b2e8cba5379cca97",
    "0x0930315cae1f9a6017c3f0c8f2314baa130e1cf13f6532bff0a8a1790cd70af918088c3db94bda214e896e1543629795",
    "0x0e6a42010cf435fb5bacc156a585e1ea3294cc81d0ceb81924d95040298380b164f702275892cedd81b62de3aba3f6b5"
};

const char *test_vec_G2_field_elems_1_1[5] = {
    "0x145a81e418d4010cc027a68f14391b30074e89e60ee7a22f87217b2f6eb0c4b94c9115b436e6fa4607e95a98de30a435",
    "0x08b852331c96ed983e497ebc6dee9b75e373d923b729194af8e72a051ea586f3538a6ebb1e80881a082fa2b24df9f566",
    "0x01897665d9cb5db16a27657760bbea7951f67ad68f8d55f7113f24ba6ddd82caef240a9bfa627972279974894701d975",
    "0x10c4df2cacf67ea3cb3108b00d4cbd0b3968031ebc8eac4b1ebcefe84d6b715fde66bef0219951ece29d1facc8a520ef",
    "0x117d9a0defc57a33ed208428cb84e54c85a6840e7648480ae428838989d25d97a0af8e3255be62b25c2a85630d2dddd8"
};

const char *test_vec_G1_mapped_points_0_x[5] = {
    "0x11a3cce7e1d90975990066b2f2643b9540fa40d6137780df4e753a8054d07580db3b7f1f03396333d4a359d1fe3766fe",
    "0x125435adce8e1cbd1c803e7123f45392dc6e326d292499c2c45c5865985fd74fe8f042ecdeeec5ecac80680d04317d80",
    "0x08834484878c217682f6d09a4b51444802fdba3d7f2df9903a0ddadb92130ebbfa807fffa0eabf257d7b48272410afff",
    "0x0cbd7f84ad2c99643fea7a7ac8f52d63d66cefa06d9a56148e58b984b3dd25e1f41ff47154543343949c64f88d48a710",
    "0x0cf97e6dbd0947857f3e578231d07b309c622ade08f2c08b32ff372bd90db19467b2563cc997d4407968d4ac80e154f8"
};

const char *test_vec_G1_mapped_points_0_y[5] = {
    "0x0eeaf6d794e479e270da10fdaf768db4c96b650a74518fc67b04b03927754bac66f3ac720404f339ecdcc028afa091b7",
    "0x0e8828948c989126595ee30e4f7c931cbd6f4570735624fd25aef2fa41d3f79cfb4b4ee7b7e55a8ce013af2a5ba20bf2",
    "0x0b318f7ecf77f45a0f038e62d7098221d2dbbca2a394164e2e3fe953dc714ac2cde412d8f2d7f0c03b259e6795a2508e",
    "0x052c00e4ed52d000d94881a5638ae9274d3efc8bc77bc0e5c650de04a000b2c334a9e80b85282a00f3148dfdface0865",
    "0x127f0cddf2613058101a5701f4cb9d0861fd6c2a1b8e0afe194fccf586a3201a53874a2761a9ab6d7220c68661a35ab3"
};

const char *test_vec_G1_mapped_points_1_x[5] = {
    "0x160003aaf1632b13396dbad518effa00fff532f604de1a7fc2082ff4cb0afa2d63b2c32da1bef2bf6c5ca62dc6b72f9c",
    "0x11def93719829ecda3b46aa8c31fc3ac9c34b428982b898369608e4f042babee6c77ab9218aad5c87ba785481eff8ae4",
    "0x158418ed6b27e2549f05531a8281b5822b31c3bf3144277fbb977f8d6e2694fedceb7011b3c2b192f23e2a44b2bd106e",
    "0x06493fb68f0d513af08be0372f849436a787e7b701ae31cb964d968021d6ba6bd7d26a38aaa5a68e8c21a6b17dc8b579",
    "0x092f1acfa62b05f95884c6791fba989bbe58044ee6355d100973bf9553ade52b47929264e6ae770fb264582d8dce512a"
};

const char *test_vec_G1_mapped_points_1_y[5] = {
    "0x0d8bb2d14e20cf9f6036152ed386d79189415b6d015a20133acb4e019139b94e9c146aaad5817f866c95d609a361735e",
    "0x0007c9cef122ccf2efd233d6eb9bfc680aa276652b0661f4f820a653cec1db7ff69899f8e52b8e92b025a12c822a6ce6",
    "0x1879074f344471fac5f839e2b4920789643c075792bec5af4282c73f7941cda5aa77b00085eb10e206171b9787c4169f",
    "0x02e98f2ccf5802b05ffaac7c20018bc0c0b2fd580216c4aa2275d2909dc0c92d0d0bdc979226adeb57a29933536b6bb4",
    "0x028e6d0169a72cfedb737be45db6c401d3adfb12c58c619c82b93a5dfcccef12290de530b0480575ddc8397cda0bbebf"
};

const char *test_vec_G2_mapped_points_0_0_x[5] = {
    "0x019ad3fc9c72425a998d7ab1ea0e646a1f6093444fc6965f1cad5a3195a7b1e099c050d57f45e3fa191cc6d75ed7458c",
    "0x12b2e525281b5f4d2276954e84ac4f42cf4e13b6ac4228624e17760faf94ce5706d53f0ca1952f1c5ef75239aeed55ad",
    "0x0f48f1ea1318ddb713697708f7327781fb39718971d72a9245b9731faaca4dbaa7cca433d6c434a820c28b18e20ea208",
    "0x09eccbc53df677f0e5814e3f86e41e146422834854a224bf5a83a50e4cc0a77bfc56718e8166ad180f53526ea9194b57",
    "0x17cadf8d04a1a170f8347d42856526a24cc466cb2ddfd506cff01191666b7f944e31244d662c904de5440516a2b09004"
};

const char *test_vec_G2_mapped_points_0_0_y[5] = {
    "0x171c88b0b0efb5eb2b88913a9e74fe111a4f68867b59db252ce5868af4d1254bfab77ebde5d61cd1a86fb2fe4a5a1c1d",
    "0x05d8a724db78e570e34100c0bc4a5fa84ad5839359b40398151f37cff5a51de945c563463c9efbdda569850ee5a53e77",
    "0x06051467c8f85da5ba2540974758f7a1e0239a5981de441fdd87680a995649c211054869c50edbac1f3a86c561ba3162",
    "0x0c3633943f91daee715277bd644fba585168a72f96ded64fc5a384cce4ec884a4c3c30f08e09cd2129335dc8f67840ec",
    "0x0d13ba91f2a8b0051cf3279ea0ee63a9f19bc9cb8bfcc7d78b3cbd8cc4fc43ba726774b28038213acf2b0095391c523e"
};

const char *test_vec_G2_mapped_points_0_1_x[5] = {
    "0x0ba10604e62bdd9eeeb4156652066167b72c8d743b050fb4c1016c31b505129374f76e03fa127d6a156213576910fef3",
    "0x02eacdc556d0bdb5d18d22f23dcb086dd106cad713777c7e6407943edbe0b3d1efe391eedf11e977fac55f9b94f2489c",
    "0x168b3d6df80069dbbedb714d41b32961ad064c227355e1ce5fac8e105de5e49d77f0c64867f3834848f152497eb76333",
    "0x0eb6186a0457d5b12d132902d4468bfeb7315d83320b6c32f1c875f344efcba979952b4aa418589cb01af712f98cc555",
    "0x17ef19497d6d9246fa94d35575c0f8d06ee02f21a284dbeaa78768cb1e25abd564e3381de87bda26acd04f41181610c5"

};

const char *test_vec_G2_mapped_points_0_1_y[5] = {
    "0x0eb22c7a543d3d376e9716a49b72e79a89c9bfe9feee8533ed931cbb5373dde1fbcd7411d8052e02693654f71e15410a",
    "0x04bbe48bfd5814648d0b9e30f0717b34015d45a861425fabc1ee06fdfce36384ae2c808185e693ae97dcde118f34de41",
    "0x134e0e8331cee8cb12f9c2d0742714ed9eee78a84d634c9a95f6a7391b37125ed48bfc6e90bf3546e99930ff67cc97bc",
    "0x119e3cf167e69eb16c1c7830e8df88856d48be12e3ff0a40791a5cd2f7221311d4bf13b1847f371f467357b3f3c0b4c7",
    "0x12c3c913ba4ed03c24f0721a81a6be7430f2971ffca8fd1729aafe496bb725807531b44b34b59b3ae5495e5a2dcbd5c8"
};

const char *test_vec_G2_mapped_points_1_0_x[5] = {
    "0x113d2b9cd4bd98aee53470b27abc658d91b47a78a51584f3d4b950677cfb8a3e99c24222c406128c91296ef6b45608be",
    "0x19f18cc5ec0c2f055e47c802acc3b0e40c337256a208001dde14b25afced146f37ea3d3ce16834c78175b3ed61f3c537",
    "0x004fd03968cd1c99a0dd84551f44c206c84dcbdb78076c5bfee24e89a92c8508b52b88b68a92258403cbe1ea2da3495f",
    "0x0eb3aabc1ddfce17ff18455fcc7167d15ce6b60ddc9eb9b59f8d40ab49420d35558686293d046fc1e42f864b7f60e381",
    "0x16ec57b7fe04c71dfe34fb5ad84dbce5a2dbbd6ee085f1d8cd17f45e8868976fc3c51ad9eeda682c7869024d24579bfd"
};

const char *test_vec_G2_mapped_points_1_0_y[5] = {
    "0x13855912321c5cb793e9d1e88f6f8d342d49c0b0dbac613ee9e17e3c0b3c97dfbb5a49cc3fb45102fdbaf65e0efe2632",
    "0x15b0dadc256a258b4c68ea43605dffa6d312eef215c19e6474b3e101d33b661dfee43b51abbf96fee68fc6043ac56a58",
    "0x1674338ea298281b636b2eb0fe593008d03171195fd6dcd4531e8a1ed1f02a72da238a17a635de307d7d24aa2d969a47",
    "0x198bdfb19d7441ebcca61e8ff774b29d17da16547d2c10c273227a635cacea3f16826322ae85717630f0867539b5ed8b",
    "0x13103f7aace1ae1420d208a537f7d3a9679c287208026e4e3439ab8cd534c12856284d95e27f5e1f33eec2ce656533b0"
};

const char *test_vec_G2_mapped_points_1_1_x[5] = {
    "0x0fd3def0b7574a1d801be44fde617162aa2e89da47f464317d9bb5abc3a7071763ce74180883ad7ad9a723a9afafcdca",
    "0x05e47c1781286e61c7ade887512bd9c2cb9f640d3be9cf87ea0bad24bd0ebfe946497b48a581ab6c7d4ca74b5147287f",
    "0x0dc7fa13fff6b12558419e0a1e94bfc3cfaf67238009991c5f24ee94b632c3d09e27eca329989aee348a67b50d5e236c",
    "0x0aaf1dee3adf3ed4c80e481c09b57ea4c705e1b8d25b897f0ceeec3990748716575f92abff22a1c8f4582aff7b872d52",
    "0x0958b2c4c2c10fcef5a6c59b9e92c4a67b0fae3e2e0f1b6b5edad9c940b8f3524ba9ebbc3f2ceb3cfe377655b3163bd7"
};

const char *test_vec_G2_mapped_points_1_1_y[5] = {
    "0x056f617902b3c0d0f78a9a8cbda43a26b65f602f8786540b9469b060db7b38417915b413ca65f875c130bebfaa59790c",
    "0x19f98db2f4a1fcdf56a9ced7b320ea9deecf57c8e59236b0dc21f6ee7229aa9705ce9ac7fe7a31c72edca0d92370c096",
    "0x169585e164c131103d85324f2d7747b23b91d66ae5d947c449c8194a347969fc6bbd967729768da485ba71868df8aed2",
    "0x0d058d9061ed27d4259848a06c96c5ca68921a5d269b078650c882cb3c2bd424a8702b7a6ee4e0ead9982baf6843e924",
    "0x0ccb594ed8bd14ca64ed9cb4e0aba221be540f25dd0d6ba15a4a4be5d67bcf35df7853b2d8dad3ba245f1ea3697f66aa"
};

const char *test_vec_G1_hash_out_x[5] = {
    "0x052926add2207b76ca4fa57a8734416c8dc95e24501772c814278700eed6d1e4e8cf62d9c09db0fac349612b759e79a1",
    "0x03567bc5ef9c690c2ab2ecdf6a96ef1c139cc0b2f284dca0a9a7943388a49a3aee664ba5379a7655d3c68900be2f6903",
    "0x11e0b079dea29a68f0383ee94fed1b940995272407e3bb916bbf268c263ddd57a6a27200a784cbc248e84f357ce82d98",
    "0x15f68eaa693b95ccb85215dc65fa81038d69629f70aeee0d0f677cf22285e7bf58d7cb86eefe8f2e9bc3f8cb84fac488",
    "0x082aabae8b7dedb0e78aeb619ad3bfd9277a2f77ba7fad20ef6aabdc6c31d19ba5a6d12283553294c1825c4b3ca2dcfe"
};

const char *test_vec_G1_hash_out_y[5] = {
    "0x08ba738453bfed09cb546dbb0783dbb3a5f1f566ed67bb6be0e8c67e2e81a4cc68ee29813bb7994998f3eae0c9c6a265",
    "0x0b9c15f3fe6e5cf4211f346271d7b01c8f3b28be689c8429c85b67af215533311f0b8dfaaa154fa6b88176c229f2885d",
    "0x03a87ae2caf14e8ee52e51fa2ed8eefe80f02457004ba4d486d6aa1f517c0889501dc7413753f9599b099ebcbbd2d709",
    "0x1807a1d50c29f430b8cafc4f8638dfeeadf51211e1602a5f184443076715f91bb90a48ba1e370edce6ae1062f5e6dd38",
    "0x05b84ae5a942248eea39e1d91030458c40153f3b654ab7872d779ad1e942856a20c438e8d99bc8abfbf74729ce1f7ac8"
};

const char *test_vec_G2_hash_out_x_a[5] = {
    "0x0141ebfbdca40eb85b87142e130ab689c673cf60f1a3e98d69335266f30d9b8d4ac44c1038e9dcdd5393faf5c41fb78a",
    "0x02c2d18e033b960562aae3cab37a27ce00d80ccd5ba4b7fe0e7a210245129dbec7780ccc7954725f4168aff2787776e6",
    "0x121982811d2491fde9ba7ed31ef9ca474f0e1501297f68c298e9f4c0028add35aea8bb83d53c08cfc007c1e005723cd0",
    "0x19a84dd7248a1066f737cc34502ee5555bd3c19f2ecdb3c7d9e24dc65d4e25e50d83f0f77105e955d78f4762d33c17da",
    "0x01a6ba2f9a11fa5598b2d8ace0fbe0a0eacb65deceb476fbbcb64fd24557c2f4b18ecfc5663e54ae16a84f5ab7f62534"
};

const char *test_vec_G2_hash_out_x_b[5] = {
    "0x05cb8437535e20ecffaef7752baddf98034139c38452458baeefab379ba13dff5bf5dd71b72418717047f5b0f37da03d",
    "0x139cddbccdc5e91b9623efd38c49f81a6f83f175e80b06fc374de9eb4b41dfe4ca3a230ed250fbe3a2acf73a41177fd8",
    "0x190d119345b94fbd15497bcba94ecf7db2cbfd1e1fe7da034d26cbba169fb3968288b3fafb265f9ebd380512a71c3f2c",
    "0x0934aba516a52d8ae479939a91998299c76d39cc0c035cd18813bec433f587e2d7a4fef038260eef0cef4d02aae3eb91",
    "0x11fca2ff525572795a801eed17eb12785887c7b63fb77a42be46ce4a34131d71f7a73e95fee3f812aea3de78b4d01569"
};

const char *test_vec_G2_hash_out_y_a[5] = {
    "0x0503921d7f6a12805e72940b963c0cf3471c7b2a524950ca195d11062ee75ec076daf2d4bc358c4b190c0c98064fdd92",
    "0x1787327b68159716a37440985269cf584bcb1e621d3a7202be6ea05c4cfe244aeb197642555a0645fb87bf7466b2ba48",
    "0x05571a0f8d3c08d094576981f4a3b8eda0a8e771fcdcc8ecceaf1356a6acf17574518acb506e435b639353c2e14827c8",
    "0x14f81cd421617428bc3b9fe25afbb751d934a00493524bc4e065635b0555084dd54679df1536101b2c979c0152d09192",
    "0x0b6798718c8aed24bc19cb27f866f1c9effcdbf92397ad6448b5c9db90d2b9da6cbabf48adc1adf59a1a28344e79d57e"
};

const char *test_vec_G2_hash_out_y_b[5] = {
    "0x12424ac32561493f3fe3c260708a12b7c620e7be00099a974e259ddc7d1f6395c3c811cdd19f1e8dbf3e9ecfdcbab8d6",
    "0x00aa65dae3c8d732d10ecd2c50f8a1baf3001578f71c694e03866e9f3d49ac1e1ce70dd94a733534f106d4cec0eddd16",
    "0x0bb5e7572275c567462d91807de765611490205a941a5a6af3b1691bfe596c31225d3aabdf15faff860cb4ef17c7c3be",
    "0x09bcccfa036b4847c9950780733633f13619994394c23ff0b32fa6b795844f4a0673e20282d07bc69641cee04f5e5662",
    "0x03a47f8e6d1763ba0cad63d6114c0accbef65707825a511b251a660a9b3994249ae4e63fac38b23da0c398689ee2ab52"
};

bool test_message_expand_xmd(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    uint32_t len_in_bytes;
    octet_string *bytes;
    char *out;

    len_in_bytes = test_vec_len_in_bytes[test_num];

    octet_string_alloc(&bytes, len_in_bytes);

    expand_message_xmd(bytes, test_vec_msg[test_num],
                       EXPAND_MESSAGE_XMD_TEST_DST, len_in_bytes);

    out = octet_string_to_str(bytes);
    ret = strcmp(out, test_vec_uniform_bytes[test_num]) == 0;

    printf("Actual:\n\t%s\n", out);
    printf("Expected:\n\t%s\n", test_vec_uniform_bytes[test_num]);

    octet_string_free(bytes);
    free(out);

    return ret;
}

bool test_hash_to_field_fp(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    mpz_t *elems;
    mpz_t expect;

    mpz_init(expect);

    elems = hash_to_field_fp(test_vec_msg[test_num], G1_TEST_DST, 2);

    mpz_set_str(expect, test_vec_G1_field_elems_0[test_num], 0);
    ret = (mpz_cmp(expect, elems[0]) == 0);
    printf("u[0]\nActual:\n\t%s\n", mpz_get_str(NULL, 16, elems[0]));
    printf("Expected:\n\t%s\n", mpz_get_str(NULL, 16, expect));

    mpz_set_str(expect, test_vec_G1_field_elems_1[test_num], 0);
    ret = ret && (mpz_cmp(expect, elems[1]) == 0);
    printf("u[1]\nActual:\n\t%s\n", mpz_get_str(NULL, 16, elems[1]));
    printf("Expected:\n\t%s\n", mpz_get_str(NULL, 16, expect));

    mpz_clear(expect);

    return ret;
}

bool test_hash_to_field_fp2(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    fp2_elem **elems;
    fp2_elem expect;

    elems = hash_to_field_fp2(test_vec_msg[test_num], G2_TEST_DST, 2);

    fp2_elem_from_str(&expect,
            test_vec_G2_field_elems_0_0[test_num],
            test_vec_G2_field_elems_0_1[test_num]);

    ret = fp2_equal(&expect, elems[0]);
    printf("u[0]\nActual:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, elems[0]->a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, elems[0]->b));
    printf("Expected:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, expect.a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, expect.b));
    fp2_elem_free(&expect);

    fp2_elem_from_str(&expect,
            test_vec_G2_field_elems_1_0[test_num],
            test_vec_G2_field_elems_1_1[test_num]);

    ret = ret && fp2_equal(&expect, elems[1]);
    printf("u[1]\nActual:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, elems[1]->a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, elems[1]->b));
    printf("Expected:\n");
    printf("\ta: %s\n", mpz_get_str(NULL, 16, expect.a));
    printf("\tb: %s\n", mpz_get_str(NULL, 16, expect.b));
    fp2_elem_free(&expect);

    fp2_elem_free(elems[0]);
    fp2_elem_free(elems[1]);
    free(elems);

    return ret;
}

bool test_map_to_curve_G1(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    mpz_t u;
    G1_elem_affine actual, expect;

    mpz_init(u);
    G1_identity_init_affine(&actual);

    mpz_set_str(u, test_vec_G1_field_elems_0[test_num], 0);
    G1_elem_affine_from_str(&expect,
                           test_vec_G1_mapped_points_0_x[test_num],
                           test_vec_G1_mapped_points_0_y[test_num]);

    map_to_curve_G1(&actual, u);

    ret = G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    G1_elem_free_affine(&expect);

    mpz_set_str(u, test_vec_G1_field_elems_1[test_num], 0);
    G1_elem_affine_from_str(&expect,
                           test_vec_G1_mapped_points_1_x[test_num],
                           test_vec_G1_mapped_points_1_y[test_num]);

    map_to_curve_G1(&actual, u);

    ret = ret && G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));
    G1_elem_free_affine(&expect);

    G1_elem_free_affine(&actual);
    mpz_clear(u);

    return ret;
}

bool test_map_to_curve_G2(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    fp2_elem u;
    G2_elem_affine actual, expect;

    fp2_elem_init(&u);
    G2_identity_init_affine(&actual);

    fp2_elem_set_str(&u,
            test_vec_G2_field_elems_0_0[test_num],
            test_vec_G2_field_elems_0_1[test_num]);
    G2_elem_affine_from_str(&expect,
            test_vec_G2_mapped_points_0_0_x[test_num],
            test_vec_G2_mapped_points_0_0_y[test_num],
            test_vec_G2_mapped_points_0_1_x[test_num],
            test_vec_G2_mapped_points_0_1_y[test_num]);

    map_to_curve_G2(&actual, &u);

    ret = G2_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));
    G2_elem_free_affine(&expect);

    fp2_elem_set_str(&u,
            test_vec_G2_field_elems_1_0[test_num],
            test_vec_G2_field_elems_1_1[test_num]);
    G2_elem_affine_from_str(&expect,
            test_vec_G2_mapped_points_1_0_x[test_num],
            test_vec_G2_mapped_points_1_0_y[test_num],
            test_vec_G2_mapped_points_1_1_x[test_num],
            test_vec_G2_mapped_points_1_1_y[test_num]);

    map_to_curve_G2(&actual, &u);

    ret = ret && G2_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));
    G2_elem_free_affine(&expect);

    G2_elem_free_affine(&actual);
    fp2_elem_free(&u);

    return ret;
}

bool test_hash_to_G1(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    G1_elem_affine actual, expect;

    G1_identity_init_affine(&actual);

    G1_elem_affine_from_str(&expect,
                           test_vec_G1_hash_out_x[test_num],
                           test_vec_G1_hash_out_y[test_num]);

    BLS12_381_hash_to_G1(&actual, test_vec_msg[test_num], G1_TEST_DST);

    ret = G1_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, actual.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, actual.y));
    printf("Expected:\n");
    printf("\tx: %s\n", mpz_get_str(NULL, 16, expect.x));
    printf("\ty: %s\n", mpz_get_str(NULL, 16, expect.y));

    G1_elem_free_affine(&expect);
    G1_elem_free_affine(&actual);

    return ret;
}

bool test_hash_to_G2(int test_num)
{
    printf("Running %s #%d\n", __func__, test_num);

    bool ret;
    G2_elem_affine actual, expect;

    G2_identity_init_affine(&actual);

    G2_elem_affine_from_str(&expect,
                           test_vec_G2_hash_out_x_a[test_num],
                           test_vec_G2_hash_out_x_b[test_num],
                           test_vec_G2_hash_out_y_a[test_num],
                           test_vec_G2_hash_out_y_b[test_num]);

    BLS12_381_hash_to_G2(&actual, test_vec_msg[test_num], G2_TEST_DST);

    ret = G2_equiv_affine(&actual, &expect);
    printf("Actual:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, actual.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, actual.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, actual.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, actual.y->b));
    printf("Expected:\n");
    printf("\tx0: %s\n", mpz_get_str(NULL, 16, expect.x->a));
    printf("\tx1: %s\n", mpz_get_str(NULL, 16, expect.x->b));
    printf("\ty0: %s\n", mpz_get_str(NULL, 16, expect.y->a));
    printf("\ty1: %s\n", mpz_get_str(NULL, 16, expect.y->b));

    G2_elem_free_affine(&expect);
    G2_elem_free_affine(&actual);

    return ret;
}

int main()
{
    TEST_MAIN_INIT;

    fp_params_init();

    for (int i = 0; i < 10; i++) {
        TEST_RUN(test_message_expand_xmd(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_field_fp(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_field_fp2(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_map_to_curve_G1(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_map_to_curve_G2(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_G1(i));
    }

    for (int i = 0; i < 5; i++) {
        TEST_RUN(test_hash_to_G2(i));
    }

    fp_params_free();

    TEST_MAIN_RETURN;
}
