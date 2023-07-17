#include <gtest/gtest.h>

#include "angle.h"

TEST(Angle, ToRadians)
{
	EXPECT_TRUE(std::abs(math::toRadians(180) - math::PI) < math::EPS);
}
TEST(Angle, ToDegrees)
{
	EXPECT_TRUE(std::abs(math::toDegrees(math::PI) - 180) < math::EPS);
}